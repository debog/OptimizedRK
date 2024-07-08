function [contour_matrix] = plotStabilityRegion(ifig,p,lambda,bounds,ls,lw)
% function [contour_matrix] = plotStabilityRegion(ifig,p,lambda,bounds,ls,lw)
%
% plot the absolute stability region of a one-step method,
% given the stability function
%
% Inputs:
%       * ifig: figure index
%       * p: coefficients of the stability function
%       * lambda: eigenvalues to plot along with stability function
%       * bounds: bounds for region to compute and plot (default [-9 1 -5 5])
%       * ls:   line style (default '-r')
%       * lw:   line width (default 2)

if nargin<6 lw=4; end
if nargin<5 ls='-r'; end
if nargin<4 bounds=[-9 1 -5 5]; end

figure(ifig);
q=[1];

finished=0;
while ~finished
    dx=(bounds(2)-bounds(1))/500.;
    xa=bounds(1):dx:bounds(2); ya=bounds(3):dx:bounds(4);

    X=repmat(xa,length(ya),1); Y=repmat(ya',1,length(xa));
    XY=(X+1i*Y);
    S=ones(size(XY)); T=ones(size(XY));
    m=length(p)-1; n=length(q)-1;

    %Evaluate numerator
    XYP=XY;
    for j=1:m
      S=S+p(j+1)*XYP; XYP=XYP.*XY;
    end
    %Evaluate denominator
    XYP=XY;
    for j=1:n
      T=T+q(j+1)*XYP; XYP=XYP.*XY;
    end

    %Compute $R(z)=|S(z)/T(z)|$
    R=abs(S./T);

    % Now check whether the whole stability region is contained in the plot.
    % If not, make it bigger.  Note that this will loop infinitely for
    % some implicit methods!
    finished=1;
    if any(R(1,:)<1)
        finished=0;
        bounds(3:4)=bounds(3:4)*1.3;
    end
    if any(R(:,1)<1)
        finished=0;
        bounds(1)=bounds(1)*1.3;
    end
    % Safety valve
    if any(bounds>1000)
        finished=1;
    end
end

%Plot the absolute stability boundary ($R=1$)
contour_matrix = contour(xa,ya,R,[1 1],ls,'LineWidth',lw);
title('Absolute stability region and scaled eigenvalues');
hold on;
axis equal
v=axis;
plot([0 0],v(3:4),'--k', 'HandleVisibility', 'off');
hold on;
plot(v(1:2),[0 0],'--k', 'HandleVisibility', 'off');
hold on;
plot(real(lambda),imag(lambda),'bo');
hold off;
