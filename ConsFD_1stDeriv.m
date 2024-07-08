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
% Boundary type
boundary = 'periodic';
% Discretization method ('1','weno5','crweno5')
sp_method = 'weno5';
% Set parameters for RK time integrator
stages = 5;
order = 2;
class = 'erk'; % explicit RK
obj_crit = 'acc'; % truncation error
num_proc = 4;
if (obj_crit == 'acc')
    opt_alg = 'interior-point';
else
    opt_alg = 'sqp';
end

%% Construct the discretization matrix
fprintf('Computing discretization matrix.\n');
InterpMatrix = GetInterpOperator(N,sp_method,boundary);
FDMatrix = GetFDOperator(N);
DiscretMatrix = -FDMatrix*InterpMatrix;

%% Compute and plot the spectrum of the RHS Jacobian matrix
figure(ifig);
fprintf('Computing spectrum.\n');
lambda = eig(DiscretMatrix);
figure(ifig);
plot(real(lambda),imag(lambda),'bo');
title('Eigenvalues of the discretization matrix');
axis equal;
grid on;
ifig = ifig + 1;

%% Compute optimized stability polynomial
cvx_clear;
tol = 1.e-2;
[cfl_max, poly_coeff] = opt_poly_bisect(lambda, stages, order, 'chebyshev');

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
    plotStabilityRegionRK(ifig,rk,cfl_max*lambda);
    ifig = ifig + 1;
else
    fprintf('Unable to find RK method.\n');
    plotStabilityRegion(ifig,poly_coeff,cfl_max*lambda);
    ifig = ifig + 1;
end
