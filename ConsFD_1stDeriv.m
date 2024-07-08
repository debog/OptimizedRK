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
% Physical model('linadv1d','euler1d')
model = 'linadv1d';
% Boundary type
boundary = 'periodic';
% Discretization method ('1','weno5','crweno5')
sp_method = 'weno5';
% Flow variables (if using 'euler1d')
gamma = 1.4;
rho_inf = 1.0;
drho = 0.1;
u_inf = 0.1;
p_inf = rho_inf/gamma;
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

%% Initialize grid
if (strcmp(strtrim(boundary),'periodic'))
    dx = 1.0/N;
else
    dx = 1.0/(N-1);
end
x = 0:dx:(N-1)*dx;

%% Initialize solution
if (strcmp(model,'linadv1d'))
    nvar = 1;
    U = sin(x);
    GetFluxJacobian = @GetFluxJacobian_LinAdv1D;
elseif (strcmp(model,'euler1d'))
    nvar = 3;
    rho = rho_inf+drho*sin(x);
    u = u_inf*ones(1,size(x,2));
    p = p_inf*ones(1,size(x,2));
    e = p/(gamma-1) + 0.5 * rho .* u.^2;
    U = [rho;rho.*u;e];
    GetFluxJacobian = @GetFluxJacobian_Euler1D;
else
    fprintf('Error: unknown physical model %s\n',model);
    return;
end
ndof = N * nvar;

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
