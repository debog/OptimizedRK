function lambda = Get3rdDerivSpectrum( N,method )
%% GET3RDDERIVSPECTRUM Returns the spectrum of the finite-difference discretization
%%                     of the 3rd derivative

if (strcmp(strtrim(method),'tdccs8'))
    lambda = TDCCS8(N);
elseif (strcmp(strtrim(method),'tdcncs8'))
    lambda = TDCNCS8(N);
else
    lambda = zeros(N+1,N);
end

%% Done
end

function lambda = TDCCS8( N )
%% TDCCS8 Returns the matrix corresponding to the TDCCS8 scheme
    theta = linspace(0,2*pi,N);
    coeff_a=58021/14120; coeff_b = -109007/28240; coeff_c=1029/28240; coeff_alpha = -1261/3530 ;
    denom = (1+2*coeff_alpha*cos(theta)) ;
    num = -(( 2*coeff_a*(8.*sin(theta/2)-4.*sin(theta)) ) + ( 2*coeff_b/5*(12.*sin(theta)-8.*sin(3*theta/2))) + ( 2*coeff_c/35*(20.*sin(theta)-8.*sin(5*theta/2))));
    lambda = ((0.*theta)+1i*(num./denom)).';
%% Done
end


function lambda = TDCNCS8( N )
%% TDCNCS8 Returns the matrix corresponding to the TDCNCS8 scheme
    theta = linspace(0,2*pi,N);
    coeff_a=2367/1180; coeff_b = -167/1180; coeff_c=1/236; coeff_alpha = 205/472 ;
    denom = (1+(2*coeff_alpha*cos(theta))) ;
    num = -(( coeff_a*(2.*sin(theta)-1.*sin(2*theta)) ) + ( coeff_b/4*(3.*sin(theta)-1.*sin(3*theta)) ) + ( coeff_c/10*(4.*sin(theta)-1.*sin(4*theta)) )  );
    lambda = ((0.*theta)+1i*(num./denom)).';
%% Done
end

