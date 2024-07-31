clc;
clear all;
close all;
ifig = 1;

% adding paths
fprintf('Adding current directory to path.\n');
addpath(genpath('./'));

fprintf('Initializing...\n');
%% Set some parameters 
% Grid size
N = 100;
% Discretization method ('tdccs8','tdcncs8')
sp_method = 'tdcncs8'; 
% Set parameters for RK time integrator
stages = 4;
order = 3;
class = 'erk'; % explicit RK
 obj_crit = 'acc'; % minimize leading trunction error coefficient
num_proc = 4;   % num_proc = 4;
if (obj_crit == 'acc')
    opt_alg = 'interior-point';
else
    opt_alg = 'sqp';
end

%% Plot the spectrum of the RHS Jacobian matrix
figure(ifig);
fprintf('Computing spectrum.\n');
lambda = Get3rdDerivSpectrum(N,sp_method);

%% Compute optimized stability polynomial
fprintf('Computing optimized stability polynomial.\n');
cvx_clear;
tol = 1.e-2;
[cfl_max, poly_coeff] = opt_poly_bisect(lambda, stages, order, 'monomial');

%% Compute ERK method with this stability polynomial
num_coeffs = stages - order;
p_ind = zeros(num_coeffs,1);
p_val = zeros(num_coeffs,1);
for i = 1:num_coeffs
    p_ind(i) = i+order;
    p_val(i) = poly_coeff(i+order+1);
end
rk = rk_opt(  stages, ...
              order, ...
              class, ...
              obj_crit, ...
              'startvec','smart', ...
              'solveorderconditions',1, ...
              'algorithm',opt_alg, ...
              'poly_coeff_ind',p_ind, ...
              'poly_coeff_val',p_val, ...
              'np', num_proc );
fprintf('Found stability polynomial with CFL = %f\n',cfl_max);
if (isstruct(rk))
    fprintf('Found optimized RK method of order %d with %d stages:\n', order, stages);
    fprintf('Butcher table:\n');
    fprintf('  A = \n');
    for i = 1:stages
        fprintf('    ');
        for j = 1:stages
            fprintf('%1.4f ', rk.A(i,j));
        end
        fprintf('\n');
    end
    fprintf('  b = \n');
    fprintf('    ');
    for j = 1:stages
        fprintf('%1.4f ', rk.b(j));
    end
    fprintf('\n');
    % fprintf('  r = %f\n',rk.r);
    plotStabilityRegionRK(ifig,rk,cfl_max*lambda);
    ifig = ifig + 1;
else
    fprintf('Unable to find RK method.\n');
    plotStabilityRegion(ifig,poly_coeff,cfl_max*lambda);
    ifig = ifig + 1;
end
