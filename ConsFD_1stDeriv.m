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
% Set order and stages for RK time integrator
stages = 5;
order = 1;

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
DiscretMatrix = -FDMatrix*InterpMatrix/dx;

%% Construct the RHS Jacobian matrix
fprintf('Computing Jacobian.\n');
RHSJac = zeros(ndof,ndof);
for i = 1:N
    for j = 1:N
        Amat = feval(GetFluxJacobian,U(:,j),gamma);
        for k = 1:nvar
            for l = 1:nvar
                RHSJac(nvar*(i-1)+k,nvar*(j-1)+l) = DiscretMatrix(i,j) * Amat(k,l);
            end
        end
    end
end

%% Compute and plot the spectrum of the RHS Jacobian matrix
figure(ifig);
fprintf('Computing spectrum.\n');
lambda = eig(RHSJac);
figure(ifig);
plot(real(lambda),imag(lambda),'bo');
title('Eigenvalues of the RHS Jacobian');
axis equal;
grid on;
ifig = ifig + 1;

%% Compute optimized stability polynomial
cvx_clear;
tol = 1.e-2;
[dt_max, poly_coeff] = opt_poly_bisect(lambda, stages, order, 'chebyshev');
fprintf('Found stability polynomial with dt_max= %1.3e, CFL = %f\n',dt_max, dt_max/dx);
plotStabilityRegion(ifig,poly_coeff,dt_max*lambda);
ifig = ifig + 1;

%% Compute ERK method with this stability polynomial
num_coeffs = stages - order;
p_ind = zeros(num_coeffs,1);
p_val = zeros(num_coeffs,1);
for i = 1:num_coeffs
    p_ind(i) = i+order;
    p_val(i) = poly_coeff(i+order+1);
end
rk = rk_opt(stages+1,order,'erk','acc','startvec','smart','solveorderconditions',1,'algorithm','interior-point','poly_coeff_ind',p_ind,'poly_coeff_val',p_val,'np',4);
if (isstruct(rk))
    fprintf('Found optimized RK method of order %d with %d stages:\n', order, stages+1);
    plotStabilityRegionRK(ifig,rk,dt_max*lambda);
    ifig = ifig + 1;
end
