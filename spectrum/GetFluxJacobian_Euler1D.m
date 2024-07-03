function A = GetFluxJacobian_Euler1D( U,gamma )
%GETFLUXJACOBIAN_EULER1D Flux Jacobian for 1D Euler equations

rho = U(1);
u = U(2)/rho;
e = U(3);
p = (gamma-1) * (e - 0.5 * rho*u*u);
a = sqrt(gamma*p/rho);

D = diag([0,abs(u+a),abs(u-a)]);

R = zeros(3,3);
L = zeros(3,3);

R(1,1) =  1;
R(1,2) =  rho/(2*a);
R(1,3) = -rho/(2*a);

R(2,1) =  u;
R(2,2) =  rho/(2*a) * (u+a);
R(2,3) = -rho/(2*a) * (u-a);

R(3,1) =  u^2/2;
R(3,2) =  rho/(2*a) * (u^2/2 + a^2/(gamma-1) + a*u);
R(3,3) = -rho/(2*a) * (u^2/2 + a^2/(gamma-1) - a*u);

L(1,1) =  (gamma-1)/(rho*a) * ((rho/a)*(-u^2/2 + a^2/(gamma-1)));
L(1,2) =  (gamma-1)/(rho*a) * (rho*u/a);
L(1,3) =  (gamma-1)/(rho*a) * (-rho/a);

L(2,1) =  (gamma-1)/(rho*a) * (u^2/2 - a*u/(gamma-1));
L(2,2) =  (gamma-1)/(rho*a) * (-u+a/(gamma-1)); 
L(2,3) =  (gamma-1)/(rho*a); 

L(3,1) =  (gamma-1)/(rho*a) * (-u^2/2 - a*u/(gamma-1)); 
L(3,2) =  (gamma-1)/(rho*a) * (u+a/(gamma-1)); 
L(3,3) = -(gamma-1)/(rho*a);

A = R*D*L;

end

