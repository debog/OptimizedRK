function [contour_matrix] = plotStabilityRegionRK(ifig,rk,lambda,plotbounds,ls,lw)
% function [contour_matrix] = plotStabilityRegionRK(rk,plotbounds,ls,lw)
%
% Plots the absolute stability region
% of a Runge-Kutta method, given the Butcher array

% Inputs:
%       * ifig: figure index
%       * rk: a Runge-Kutta method
%       * lambda: eigenvalues to plot along with stability function

% Remaining inputs are optional:
%       * plotbounds: bounds for region to compute and plot (default [-9 1 -5 5])
%       * ls:   line style (default '-r')
%       * lw:   line width (default 2)

if nargin<6 lw=4; end
if nargin<5 ls='-r'; end
if nargin<4 plotbounds=[-9 1 -5 5]; end

[p,q]=rk_stabfun(rk);

figure(ifig);
contour_matrix = plotstabreg_func(p,q,plotbounds,ls,lw);
title('Absolute stability region and scaled eigenvalues');
hold on;
plot(real(lambda),imag(lambda),'bo');
axis equal
hold off;
