function [k_star, Delta0_star] = optimal_scurve_k_asym_brute_force(Pgmin, Pgmax, Pg0, eta)
%% Use brute force evaluation to compute the best values for the S-curve parameters
% This is inefficient but much more straightforward, and we don't need
% efficiency here. Bisection or a Newton method would be much more
% computationally efficient, but I don't think this will be the bottleneck
% for this code overall.

Delta = -2:0.01:2;
k = 0.001:0.001:10;
k = k*eta / (Pgmax - Pgmin);
err = nan(length(k),1);

% figure;

for i=1:length(k)
    
    % Compute the S-curve for this value of k
    Delta0 = (1/k(i))*log( (Pgmax-Pg0)/(Pg0-Pgmin) );
    Pg = Pgmin + (Pgmax-Pgmin) ./ (1+exp(-k(i)*(Delta-Delta0)));
    
    % Compute the error with respect to the true piecewise linear curve
    Delta1 = (Pgmin - Pg0)/eta;
    Delta2 = (Pgmax - Pg0)/eta;

    err(i) = sum(abs(Pg(Delta > Delta1 & Delta < Delta2) - (Pg0 + eta*Delta(Delta > Delta1 & Delta < Delta2)))) + sum(abs(Pg(Delta < Delta1) - Pgmin)) + sum(abs(Pg(Delta > Delta2) - Pgmax));
    % fprintf('%i of %i\n',i,length(k));
    
    
    % clf
    % plot(Delta,Pg,'b-','LineWidth',3);
    % hold on
    % plot(Delta,Pgmin*ones(length(Delta),1),'r--');
    % plot(Delta,Pgmax*ones(length(Delta),1),'r--');
    % plot(Delta,eta*Delta + Pg0,'r--');

end

[~,idx] =  min(err);
k_star = k(idx);
Delta0_star = (1/k_star)*log( (Pgmax-Pg0)/(Pg0-Pgmin) );
