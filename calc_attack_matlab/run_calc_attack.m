%% Calculate the worst-case line flows and voltage magnitudes for a given set of compromised generators
% January 15, 2026

% Notes:
% - Very strange: the reactive power demands in the Hawaii test case are
%   apparently all zero? That doesn't make sense. Was this test case
%   designed for DC power flow only?

clear
clc
close all

define_constants;

verbose = false; % Display solver output?

%% Load test case
% mpc=loadcase('case30');
mpc = loadcase('hawaii40_with_solar');

% Editing the solar capacity
% mpc.gen(46:49, 2:5) = mpc.gen(46:49, 2:5) * 2;
% mpc.gen(46:49, 9) = mpc.gen(46:49, 9) * 2;
% Remove generators that are turned off
gen_off = mpc.gen(:,GEN_STATUS)==0;
mpc.gen(gen_off,:) = [];
mpc.gencost(gen_off,:) = [];

nbus = size(mpc.bus,1);
ngen = size(mpc.gen,1);
nbranch = size(mpc.branch,1);

Sbase = mpc.baseMVA;

%% Choose compromised generators
c = false(ngen,1);
c(end-3:end) = true;

% Set the lower generation bounds of the compromised generators to zero to 
% model the ability to shut off the generators.
mpc.gen(logical(c),PMIN) = 0;

%% Make test case Qd reasonable
% Set the loads to have reactive power demand based on a specified power
% factor. It is weird that they all have Qd = 0.
pf = 0.96; % Assumed load power factor (lagging)
mpc.bus(:,QD) = mpc.bus(:,PD) .* sqrt((1/pf^2) - 1);

%% Run OPF to get nominal values (prior to attack) and store the relevant values
res0 = runopf(mpc);

mpc.bus(:,VM) = res0.bus(:,VM);
mpc.bus(:,VA) = res0.bus(:,VA);
mpc.gen(:,PG) = res0.gen(:,PG);
mpc.gen(:,QG) = res0.gen(:,QG);
mpc.gen(:,VG) = mpc.bus(mpc.gen(:,GEN_BUS),VM);

mpc.branch(:,PF) = res0.branch(:,PF);
mpc.branch(:,QF) = res0.branch(:,QF);

mpc.branch(:,PT) = res0.branch(:,PT);
mpc.branch(:,QT) = res0.branch(:,QT);

%% Set S-curve parameters
pvpq_k = 100; % Slope parameter for the smoothed PV/PQ switching characteristic

% Set generator participation factors based on nominal outputs. These will
% be adjusted internally to the function "calc_attack" to turn off AGC for
% generators that near to their limits as these generators cannot
% participate in AGC.
alpha = mpc.gen(:,PG) ./ sum(mpc.gen(:,PG));

%% Find worst-case outcomes for the given set of compromised generators

% Placeholder variables to store results
Ift = nan(nbranch,1);
Itf = nan(nbranch,1);
Vlow = nan(nbus,1);
Vhigh = nan(nbus,1);
% 
% % Maximize flows in f->t direction
% for i=1:nbranch
%     wft = zeros(nbranch,1);
%     wtf = zeros(nbranch,1);
%     wv = zeros(nbus,1);
% 
%     wft(i) = 1;
%     fprintf('Computing Worst-Case Flow (from -> to direction) for line %i of %i\n',i,nbranch);
%     res = calc_attack(mpc, wft, wtf, wv, c, alpha, pvpq_k,verbose);
%     if res.success
%         Ift(i) = sqrt(res.f) * 100 / (mpc.branch(i, 6)/mpc.baseMVA) ;
%     end
% end
% 
% % Maximize flows in t->f direction
% for i=1:nbranch
%     wft = zeros(nbranch,1);
%     wtf = zeros(nbranch,1);
%     wv = zeros(nbus,1);
%     wtf(i) = 1;
%     fprintf('Computing Worst-Case Flow (to -> from direction) for line %i of %i\n',i,nbranch);
%     res = calc_attack(mpc, wft, wtf, wv, c, alpha, pvpq_k,verbose);
%     if res.success
%         Itf(i) = sqrt(res.f) * 100 / (mpc.branch(i, 6) / mpc.baseMVA) ;
%     end
% end
% max_violations = max(Ift, Itf);

% Minimize voltages
for i=1:nbus
    wft = zeros(nbranch,1);
    wtf = zeros(nbranch,1);
    wv = zeros(nbus,1);

    wv(i) = -1;
    fprintf('Computing Worst-Case Lower Voltage for bus %i of %i\n',i,nbus);
    res = calc_attack(mpc, wft, wtf, wv, c, alpha, pvpq_k,verbose);
    if res.success
        Vlow(i) = -res.f;
        disp(Vlow)
    end
end
% 
% % Maximize voltages
% for i=1:nbus
%     wft = zeros(nbranch,1);
%     wtf = zeros(nbranch,1);
%     wv = zeros(nbus,1);
% 
%     wv(i) = 1;
%     fprintf('Computing Worst-Case Upper Voltage for bus %i of %i\n',i,nbus);
%     res = calc_attack(mpc, wft, wtf, wv, c, alpha, pvpq_k,verbose);
%     if res.success
%         Vhigh(i) = res.f;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the results - Line violations
figure
axL = axes;
hold(axL, 'on')

% Histogram (medium grey bins), do not add to legend
hL = histogram(axL, max_violations, 10, ...
    'FaceColor', [0.5 0.5 0.5], 'EdgeColor','k', 'DisplayName','', 'HandleVisibility','off');

% Determine x-range for shading
xl = xlim(axL);
yl = ylim(axL);                         % used for patch height

% Normal (blue) region: from xl(1) to 100
pNormalL = patch(axL, [xl(1) 100 100 xl(1)], [0 0 yl(2) yl(2)], ...
    'b', 'FaceAlpha', 0.12, 'EdgeColor','none', 'DisplayName','Normal operation');

% Violation (red) region: from 100 to xl(2)
pViolL = patch(axL, [100 xl(2) xl(2) 100], [0 0 yl(2) yl(2)], ...
    'r', 'FaceAlpha', 0.12, 'EdgeColor','none', 'DisplayName','Constraint violation');

% Prominent threshold line (red dashed)
xline(axL, 100, '--r', 'LineWidth', 2.0);

% Legend: only show the two region patches
legend(axL, [pNormalL pViolL], {'Normal operation','Constraint violation'}, 'Location','northeast');
title(axL, 'Worst-Case Line Loading','FontSize',14)
xlabel(axL, 'Percentage loading (%)','FontSize',14)
ylabel(axL, 'Number of branches','FontSize',14)
grid(axL, 'on')
hold(axL, 'off')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the results - Voltage violations
% Split worst-case voltages into two stacked plots with shaded regions
figure
t = tiledlayout(2,1, 'Padding','compact', 'TileSpacing','compact');

% Common grey color for histograms
histGrey = [0.5 0.5 0.5];

% Top: Overvoltage (Vhigh)
ax1 = nexttile;
hold(ax1, 'on')
h1 = histogram(ax1, Vhigh, 20, 'FaceColor', histGrey, 'EdgeColor', 'k', ...
    'DisplayName','', 'HandleVisibility','off');

% x limits and y limits for patches
x1 = xlim(ax1); y1 = ylim(ax1);

% Blue normal region between 1.00 and 1.05 (for overvoltage normal is 1.00-1.05)
pBlue1 = patch(ax1, [1.00 1.05 1.05 1.00], [0 0 y1(2) y1(2)], ...
    [0.3 0.6 1.0], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName','Normal operation');

% Red violation region: > 1.05
pRed1 = patch(ax1, [1.05 x1(2) x1(2) 1.05], [0 0 y1(2) y1(2)], ...
    [1 0.4 0.4], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName','Constraint violation');

% Prominent red dashed limit line at 1.05
xline(ax1, 1.05, '--r', 'LineWidth', 2.0);

% Nominal vertical line at 1.00 (solid blue, thinner)
xline(ax1, 1.00, '-b', 'LineWidth', 1.2);

% Legend (upper left)
legend(ax1, [pBlue1 pRed1], {'Normal operation','Constraint violation'}, 'Location','northwest')

title(ax1, 'Worst-Case Overvoltage','FontSize',16)
xlabel(ax1, 'Voltage (p.u.)','FontSize',14)
ylabel(ax1, 'Number of buses','FontSize',14)
grid(ax1, 'on')
hold(ax1, 'off')

% Bottom: Undervoltage (Vlow)
ax2 = nexttile;
hold(ax2, 'on')
h2 = histogram(ax2, Vlow, 20, 'FaceColor', histGrey, 'EdgeColor', 'k', ...
    'DisplayName','', 'HandleVisibility','off');

x2 = xlim(ax2); y2 = ylim(ax2);

% Red violation region: < 0.95
pRed2 = patch(ax2, [x2(1) 0.95 0.95 x2(1)], [0 0 y2(2) y2(2)], ...
    [1 0.4 0.4], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName','Constraint violation');

% Blue normal region: 0.95 to 1.00
pBlue2 = patch(ax2, [0.95 1.00 1.00 0.95], [0 0 y2(2) y2(2)], ...
    [0.3 0.6 1.0], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName','Normal operation');

% Prominent red dashed line at 0.95 (min limit)
xline(ax2, 0.95, '--r', 'LineWidth', 2.0);

% Nominal vertical line at 1.00
xline(ax2, 1.00, '-b', 'LineWidth', 1.2);

% Legend (upper right)
legend(ax2, [pBlue2 pRed2], {'Normal operation','Constraint violation'}, 'Location','northeast')

title(ax2, 'Worst-Case Undervoltage','FontSize',16)
xlabel(ax2, 'Voltage (p.u.)', 'FontSize',14)
ylabel(ax2, 'Number of buses', 'FontSize', 14)
grid(ax2, 'on')
hold(ax2, 'off')

% Shared title for tiled layout
title(t, 'Worst-Case Voltage Analysis — Distributions','FontSize', 18)

% Vlow_0.96 =
% 
%     1.0726
%     1.0545
%     1.0472
%     1.0358
%     1.0639
%     1.0395
%     1.0382
%     1.0372
%     1.0596
%     1.0387
%     1.0291
%     1.0330
%     1.0349
%     1.0540
%     1.0352
%     1.0378
%     1.0282
%     1.0283
%     1.0597
%     1.0648
%     1.0393
%     1.0769
%     1.0748
%     1.0576
%     1.0784
%     1.0588
%     1.0859
%     1.0863
%     1.0763
%        NaN
%     1.0695
%     1.0683
%     1.0643
%     1.0763
%     1.0860
%     1.0614
%     1.0853