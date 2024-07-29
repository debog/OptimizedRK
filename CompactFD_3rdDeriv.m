clc;
clear all;
close all;
tic;
ifig = 1;

% adding paths
fprintf('Adding current directory to path.\n');
addpath(genpath('./'));

fprintf('Initializing...\n');
%% Set some parameters 
% Grid size
N = 100;
submethod = 1 ; % 1 for TDCCS; 2 for TDCNCS

% Set parameters for RK time integrator
stages = 4;
order = 3;
class = 'erk'; % explicit RK
% obj_crit = 'acc'; % minimize leading trunction error coefficient
obj_crit = 'ssp'; % maximize SSP coefficient
num_proc = 4;   % num_proc = 4;
if (obj_crit == 'acc')
    opt_alg = 'interior-point';
else
    opt_alg = 'sqp';
end

%% Plot the spectrum of the RHS Jacobian matrix
figure(ifig);
fprintf('Computing spectrum.\n');
switch submethod
 case 1
% Order 8 TDCCS
    theta = linspace(0,2*pi,N);
    coeff_a=58021/14120; coeff_b = -109007/28240; coeff_c=1029/28240; coeff_alpha = -1261/3530 ;
    denom = (1+2*coeff_alpha*cos(theta)) ;
    num = -(( 2*coeff_a*(8.*sin(theta/2)-4.*sin(theta)) ) + ( 2*coeff_b/5*(12.*sin(theta)-8.*sin(3*theta/2))) + ( 2*coeff_c/35*(20.*sin(theta)-8.*sin(5*theta/2))));
    lambda = ((0.*theta)+1i*(num./denom)).';
    % lambda = 1i*linspace(-147.168,0,N).'; % or use this
    disp('TDCCS-T8 ')
 case 2
% Order 8 TDCNCS
    theta = linspace(0,2*pi,N);
    coeff_a=2367/1180; coeff_b = -167/1180; coeff_c=1/236; coeff_alpha = 205/472 ;
    denom = (1+(2*coeff_alpha*cos(theta))) ;
    num = -(( coeff_a*(2.*sin(theta)-1.*sin(2*theta)) ) + ( coeff_b/4*(3.*sin(theta)-1.*sin(3*theta)) ) + ( coeff_c/10*(4.*sin(theta)-1.*sin(4*theta)) )  );
    lambda = ((0.*theta)+1i*(num./denom)).';
    disp('TDCNCS-T8 ')
 otherwise
    disp('Unknown kind of compact scheme, exiting the code...')
    return
end     

% figure(ifig);
% plot(real(lambda),imag(lambda),'bo');
% title('Eigenvalues of the discretization matrix');
% axis equal;
% grid on;
% ifig = ifig + 1;

%% Compute optimized stability polynomial
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
toc
