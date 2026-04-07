function [res] = calc_attack(mpc, wft, wtf, wv, c, alpha, pvpq_k, verbose)
%% Compute worst-case line flows or voltage magnitudes for a given set of compromised generators
%
% Inputs:
% - mpc: Matpower case, augmented with rows for the new solar generators.
% - wft: Binary vector weighting the objective function for maximizing the flows (from bus --> to bus direction)
% - wtf: Binary vector weighting the objective function for maximizing the flows (to bus --> from bus direction)
% - wv: Binary vector weighting the objective function for maximizing/minimizing the voltage magnitudes.
% - c: Binary vector of compromised generators.
% - alpha: Vector of participation factors for the generators.
% - pvpq_k: Slope parameter for the smoothed PV/PQ switching characteristic.
%
% Outputs:
% - res: Matpower format results from solving the optimization problem
%
% Notes:
% - Assumes that the Pmax for the compromised generators is equivalent to
%   Smax (inverter maximum apparent power output limit).
% - Assumes that the input mpc variable corresponds to a solved OPF problem
%   to get the generators' nominal Pg and Vstar setpoints.
% - There might be some annoyances if there are multiple generators at a
%   bus that have different nominal voltages (generators fight each other
%   to keep voltage at different nominal values). This isn't currently
%   handled, so be aware that it might be a source of bugs for weird
%   inputs.
% - The code adjusts generators at their limits to not participate in AGC.
%   They can't move in one direction so the S-curve fit is really poor,
%   leading to infeasibilities. This function finds generators near their
%   limits, set their participation factors to zero, and renormalizes the
%   participation factors. This should arguably be done outside of this
%   function, but is included here for now.

tic;

dbg = false; % Debug mode

pvpq_switching_on = true;
distslack_thresold_on = true;

mpc = loadcase(mpc);
mpc = ext2int(mpc);

sdpopts = sdpsettings('solver', 'fmincon','verbose',verbose);
sdpopts.fmincon.ScaleProblem = true;
sdpopts.fmincon.BarrierParamUpdate = 'monotone'; % 'predictor-corrector'
sdpopts.fmincon.MaxFunEvals = 1e5;
sdpopts.fmincon.MaxIter = 1000;
sdpopts.fmincon.ScaleProblem = 'obj-and-constr'; % 'none'
sdpopts.fmincon.TolCon = 1e-5;
sdpopts.fmincon.TolFun = 1e-5;
sdpopts.fmincon.TolFunVal = 1e-5;

define_constants;

V_lv = 0.6; % Lower bound on voltage magnitudes for screening out low-voltage power flow solutions.


%% Load mpc data

nbus = size(mpc.bus,1);
ngen = size(mpc.gen,1);
nbranch = size(mpc.branch,1);

Sbase = mpc.baseMVA;

Pd = mpc.bus(:,PD) / Sbase;
Qd = mpc.bus(:,QD) / Sbase;

Pgmax = mpc.gen(:,PMAX)/Sbase;
Pgmin = mpc.gen(:,PMIN)/Sbase;
Pg0 = mpc.gen(:,PG)/Sbase; % Nominal generator output
Qg0 = mpc.gen(:,QG)/Sbase; % Nominal generator output
Vstar = mpc.gen(:,VG); % Nominal generator voltage setpoint

Qgmax = mpc.gen(:,QMAX)/Sbase;
Qgmin = mpc.gen(:,QMIN)/Sbase;

% For improved numerics, bump the limits in the case where both Vg and Q
% are at at limit (within 1% of limit). Otherwise we get weird corner results.
% Note: This assumes that Qmin < 0 and Qmax > 0.
Qgmin(abs(Qg0 - Qgmin)./abs(Qgmin) <= 1e-2) = Qgmin(abs(Qg0 - Qgmin)./abs(Qgmin) <= 1e-2) * 1.05;
Qgmax(abs(Qg0 - Qgmax)./abs(Qgmax) <= 1e-2) = Qgmax(abs(Qg0 - Qgmax)./abs(Qgmax) <= 1e-2) * 1.05;

% For generators very near to their active power limits, turn off AGC model
% as they can't move in one direction). To participate in AGC, a generator 
% needs room to move both up and down.
alpha(abs(Pg0 - Pgmin)./abs(Pgmin) <= 1e-2 | abs(Pg0 - Pgmax)./abs(Pgmax) <= 1e-2 | (abs(Pg0 - Pgmin) <= 1e-2 & abs(Pgmin) <= 1e-2)) = 0;
alpha = alpha ./ sqrt(sum(alpha.^2)); % Renormalize participation factors

ref = find(mpc.bus(:,BUS_TYPE) == REF);


%% Create variables
V = sdpvar(nbus,1);
theta = sdpvar(nbus,1);
Pg = sdpvar(ngen,1);
Qg = sdpvar(ngen,1);
Delta = sdpvar(1,1);

Vg = V(mpc.gen(:,GEN_BUS));


%% Define power flows as a function of the voltage magnitudes and angles
% Adapted from my QC relaxation code.

f = mpc.branch(:,F_BUS);
t = mpc.branch(:,T_BUS);
          
Wf = V(f).^2;
Wt = V(t).^2;

Wft = V(f).*V(t).*cos(theta(f)-theta(t)) + 1i*V(f).*V(t).*sin(theta(f)-theta(t));
Wtf = conj(Wft);

% I think this was old code for getting the phase angle difference limits
% right, which aren't relevant here, so I'm commenting out.
% maxft = max(mpc.branch(:,[F_BUS T_BUS]),[],2);
% flipped = mpc.branch(:,1) == maxft(:,1);
% Wft_temp = Wft;
% Wft(flipped) = Wtf(flipped);
% Wtf(flipped) = Wft_temp(flipped);

% Line parameters
zl = mpc.branch(:,BR_R)+1i*mpc.branch(:,BR_X);
gl = real( 1 ./ zl); % Real part of line admittance
bl = imag( 1 ./ zl); % Imaginary part of line admittance
bsl= mpc.branch(:,BR_B); % Line shunt susceptance

gl = gl(:);
bl = bl(:);
yl = gl + 1i*bl;

tau = mpc.branch(:,TAP);
tau(tau == 0) = 1;
theta_sh = mpc.branch(:,SHIFT)*pi/180;
tl = tau.*exp(-1i*theta_sh(:));

% Construct flow equations
Sft = (conj(yl)-1i*bsl/2)./tau.^2 .* Wf - conj(yl./tl) .* Wft;
Stf = (conj(yl)-1i*bsl/2)         .* Wt - conj(yl)./tl .* Wtf;

p_ft = real(Sft);
p_tf = real(Stf);

q_ft = imag(Sft);
q_tf = imag(Stf);


%% Relate flows to total injections at each bus

fidx = zeros(nbus,nbranch);
tidx = zeros(nbus,nbranch);
for k=1:nbus
    % Sum power injections to get bus injections
    fidx(k,:) = mpc.branch(:,F_BUS) == k;
    tidx(k,:) = mpc.branch(:,T_BUS) == k;
end
[fidx_r,fidx_c] = find(fidx);
[tidx_r,tidx_c] = find(tidx);


rcv_f = sortrows([fidx_r fidx_c],1);
rcv_t = sortrows([tidx_r tidx_c],1);

% Adjust the columns to speed up the summation
nincident_f = sum(sparse(fidx_r, fidx_c, ones(length(fidx_r),1), nbus, nbranch),2); % number of incident lines at each bus (f_bus)
nincident_t = sum(sparse(tidx_r, tidx_c, ones(length(tidx_r),1), nbus, nbranch),2); % number of incident lines at each bus (t_bus)

max_nincident_f = max(nincident_f);
max_nincident_t = max(nincident_t);

fidx_c_adj = zeros(length(fidx_r),1);
tidx_c_adj = zeros(length(tidx_r),1);

fcount = 1;
tcount = 1;
for i=1:nbus
    fidx_c_adj(fcount:fcount+nincident_f(i)-1,1) = (1:nincident_f(i)).';
    fcount = fcount + nincident_f(i);
    
    tidx_c_adj(tcount:tcount+nincident_t(i)-1,1) = (1:nincident_t(i)).';
    tcount = tcount + nincident_t(i);
end

p_ftmat = sparse(rcv_f(:,1), fidx_c_adj, p_ft(rcv_f(:,2)), nbus, max_nincident_f);
q_ftmat = sparse(rcv_f(:,1), fidx_c_adj, q_ft(rcv_f(:,2)), nbus, max_nincident_f);

p_tfmat = sparse(rcv_t(:,1), tidx_c_adj, p_tf(rcv_t(:,2)), nbus, max_nincident_t);
q_tfmat = sparse(rcv_t(:,1), tidx_c_adj, q_tf(rcv_t(:,2)), nbus, max_nincident_t);
  
% Active and reactive power injections at each bus
P = sum(p_ftmat, 2) + sum(p_tfmat, 2) + V.^2.*mpc.bus(:,GS)/Sbase;
Q = sum(q_ftmat, 2) + sum(q_tfmat, 2) - V.^2.*mpc.bus(:,BS)/Sbase;


%% Enforce constraints

constraints = [];

% Angle reference
constraints = [constraints;
        (theta(ref) == 0):'Angle Ref';
    ];


% Ensure active power balance at each bus
Pgsum = sum(sparse(mpc.gen(:,GEN_BUS),(1:ngen).',Pg,nbus,ngen),2);
constraints = [constraints;
        (Pgsum - Pd == P):'Pg_Pinj';
    ];

% Ensure reactive power balance at each bus
Qgsum = sum(sparse(mpc.gen(:,GEN_BUS),(1:ngen).',Qg,nbus,ngen),2);
constraints = [constraints;
        (Qgsum - Qd == Q):'Qg_Qinj';
    ];



if distslack_thresold_on
    % S-curve model for non-compromised generators' active power outputs
    % (smoothed saturated participation factor model).
    for i=1:ngen
        if ~c(i)
            if alpha(i) > 0
                % Construct S-curve parameters that minimize the error with respect
                % to the piecewise linear active generation characteristic
                [k, Delta0] = optimal_scurve_k_asym_brute_force(Pgmin(i), Pgmax(i), Pg0(i), alpha(i));
                constraints = [constraints;
                        (Pg(i) == Pgmin(i) + (Pgmax(i)-Pgmin(i)) / (1+exp(-k*(Delta-Delta0)))):['Pg_Scurve gen ' int2str(i)];
                    ];
                if dbg
                    fprintf('Constructing smoothed distributed slack function for generator %i of %i\n',i,ngen);

                    figure;
                    Delta_vals = -2:0.01:2;
                    plot(Delta_vals,Pgmin(i) + (Pgmax(i)-Pgmin(i)) ./ (1+exp(-k*(Delta_vals-Delta0))),'b-','LineWidth',3);
                    hold on
                    plot(Delta_vals,Pgmin(i)*ones(length(Delta_vals),1),'r--');
                    plot(Delta_vals,Pgmax(i)*ones(length(Delta_vals),1),'r--');
                    plot(Delta_vals,alpha(i)*Delta_vals + Pg0(i),'r--');
                    plot(0,Pg0(i),'g.','MarkerSize',10);
                    title(['Smoothed DistSlack Generator ' int2str(i)]);
                end
            else
                % Generator not participating in AGC
                constraints = [constraints;
                        (Pg(i) == Pg0(i)):['Non-AGC gen ' int2str(i)];
                    ];
            end
        end
    end
else
    % Simple distributed slack bus model for non-compromised generators'
    % active power outputs (linear function, no saturation).
    for i=1:ngen
        if ~c(i)
            % Construct S-curve parameters that minimize the error with respect
            % to the piecewise linear active generation characteristic
            constraints = [constraints;
                    (Pg(i) == Pg0(i) + alpha(i)*Delta):['Simple DistSlack gen ' int2str(i)];
                ];
        end
    end
end

if pvpq_switching_on
    % S-curve model for non-compromised generators' voltage-reactive
    % characteristics (smoothed PV/PQ switching model).
    for i=1:ngen
        if ~c(i)
            constraints = [constraints;
                    % This constraint puts the S-curve through the middle of the vertical line segment
                    % (Qg(i) == Qgmin(i) + (Qgmax(i) - Qgmin(i)) / (1+exp(pvpq_k*(Vg(i)-Vstar(i))))):['PVPQ_Scurve gen ' int2str(i)];
                    % 
                    % This constraint puts the S-curve through the nominal reactive power output of the generator
                    (Qg(i) == Qgmin(i) + (Qgmax(i) - Qgmin(i)) ./ (1+exp(pvpq_k*(Vg(i)-Vstar(i)) + log( (Qgmax(i)-Qg0(i))/(Qg0(i)-Qgmin(i)) )))):['PVPQ_Scurve gen ' int2str(i)];
                ];

            if dbg
                figure;
                V_vals = 0:0.001:1.2;
                plot(V_vals,Qgmin(i) + (Qgmax(i) - Qgmin(i)) ./ (1+exp(pvpq_k*(V_vals-Vstar(i)) + log( (Qgmax(i)-Qg0(i))/(Qg0(i)-Qgmin(i)) ))),'b-','LineWidth',3);
                hold on
                plot(V_vals,Qgmin(i)*ones(length(V_vals),1),'r--');
                plot(V_vals,Qgmax(i)*ones(length(V_vals),1),'r--');
                plot(Vstar(i)*ones(2,1),[Qgmax(i); Qgmin(i)],'r--');
                plot(Vstar(i),Qg0(i),'g.','MarkerSize',10);
                title(['PVPQ Switching S-Curve: Generator ' int2str(i)]);
            end

        end
    end
else
    % PV bus generator model (constant voltage regardless of reactive power outputs).
    constraints = [constraints;
            (Vg == Vstar):'Vg_fixed';
        ];
end

% For compromised generators, enforce inverter power limits on apparent
% power as well as Pgmax, Pgmin, Qgmax, Qgmin
% Assumes that the generators' active power limit (mpc.gen(:,PMAX)) specify
% the apparent power limit for the generator.
for i=1:ngen
    if c(i)
        constraints = [constraints;
                (Pg(i)^2 + Qg(i)^2 <= Pgmax(i)^2):['Smax gen ' int2str(i)];
                (Pg(i) >= Pgmin(i)):['Pmin gen ' int2str(i)];
                % (Pg(i) <= Pgmax(i)):['Pmax gen ' int2str(i)];
                % (Qg(i) >= Qgmin(i)):['Qmin gen ' int2str(i)];
                % (Qg(i) <= Qgmax(i)):['Qmax gen ' int2str(i)];
            ];
    end
end

% Lower voltage limit in an attempt to screen out low-voltage power flow
% solutions. See footnote 1 of https://ieeexplore.ieee.org/document/8442889
constraints = [constraints;
        (V >= V_lv):'V_LV';
    ];

%% Objective function

% Define expressions for squared current flows.
Ift_sq =    V(f).^2./tau.^4 .* (gl.^2 + (bl+bsl/2).^2) ...
          + V(t).^2./tau.^2 .* (gl.^2 + bl.^2) ...
          - 2*V(f).*V(t)./tau.^3 .* ( (gl.^2 + bl.*(bl + bsl/2)).*cos(theta(f)-theta(t)-theta_sh) - gl.*bsl/2.*sin(theta(f)-theta(t)-theta_sh) );

Itf_sq =    V(t).^2 .* (gl.^2 + (bl+bsl/2).^2) ...
          + V(f).^2./tau.^2 .* (gl.^2 + bl.^2) ...
          - 2*V(f).*V(t)./tau .* ( (gl.^2 + bl.*(bl + bsl/2)).*cos(theta(f)-theta(t)-theta_sh) + gl.*bsl/2.*sin(theta(f)-theta(t)-theta_sh) );

% Construct the objective function based on the input weights.
% Negative signs convert to maximize by default.
obj = (-wft).'*Ift_sq + (-wtf).'*Itf_sq + (-wv).'*V;


%% Debug

if dbg
    % Drop in the nominal OPF solution. This should be feasible for the
    % constraints (or at least close for the smoothed PV-PQ switching).
    mpc.bus(:,VA) = mpc.bus(:,VA) - mpc.bus(ref,VA); % Set reference angle
    Vres = mpc.bus(:,VM).*exp(1i*mpc.bus(:,VA)*pi/180);
    
    assign(V,abs(Vres));
    assign(theta,angle(Vres));
    assign(Pg,mpc.gen(:,PG)/Sbase);
    assign(Qg,mpc.gen(:,QG)/Sbase);
    assign(Delta,0);

    check(constraints);
end

%% Solve

% Run solver
sdpinfo = solvesdp(constraints, obj, sdpopts);

if sdpinfo.problem == 2 || sdpinfo.problem == -3
    error(yalmiperror(sdpinfo.problem));
end
tsolve = sdpinfo.solvertime;

%% Prepare outputs

res = mpc;

res.f = double(-obj);

res.bus(:,VM) = double(V);
res.bus(:,VA) = double(theta);
res.gen(:,PG) = Pg * Sbase;
res.gen(:,QG) = Qg * Sbase;
res.gen(:,VG) = res.bus(res.gen(:,GEN_BUS),VM);

res.branch(:,PF) = double(p_ft) * Sbase;
res.branch(:,QF) = double(q_ft) * Sbase;

res.branch(:,PT) = double(p_tf) * Sbase;
res.branch(:,QT) = double(q_tf) * Sbase;

res.solinfo = sdpinfo;

res.success = sdpinfo.problem == 0;

res.tsolve = tsolve;
res.et = toc;