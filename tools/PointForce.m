
function [Ue_x,Ue_y,Un_x,Un_y,Exx_x,Exy_x,Eyy_x,Exx_y,Exy_y,Eyy_y,omega_x,omega_y] = PointForce(xs,ys,x,y,nu)

%(xs,ys) is position of source
%(x,y) are data/observation coordinates
%nu is Poisson's ratio
%Ue_x is East component due to force in x-direction
%Ue_y is East component due to force in y-direction
%Un_x is North component due to force in x-direction
%Un_y is North component due to force in y-direction
%
% Exx_j, Exy_j, Eyy_j are strain components due to force in j-direction
%
%solution for displacements in a thin elastic sheet following 
% Sandwell, D. T., & Wessel, P. (2016). Interpolation of 2?D vector data using constraints from elasticity. Geophysical Research Letters, 43(20), 10-703.


x = x-xs;
y = y-ys;
r = sqrt(x.^2 + y.^2); 

q = (3-nu)*log(r) + (1+nu)*y.^2./r.^2;
p = (3-nu)*log(r) + (1+nu)*x.^2./r.^2;
w = -(1+nu)*x.*y./r.^2;

Ue_x = q;
Ue_y = w;

Un_x = w;
Un_y = p;


%strain rates -- computed from analytical derivatives of above expressions
dq_dx = (3-nu)*x./r.^2 - 2*(1+nu)*y.^2.*x./r.^4;
dq_dy = (3-nu)*y./r.^2 + 2*(1+nu)*(-y.^3./r.^4 + y./r.^2);

dp_dx = (3-nu)*x./r.^2 + 2*(1+nu)*(-x.^3./r.^4 + x./r.^2);
dp_dy = (3-nu)*y./r.^2 - 2*(1+nu)*x.^2.*y./r.^4;

dw_dx = -(1+nu)*(-2*x.^2.*y./r.^4 + y./r.^2);
dw_dy = -(1+nu)*(-2*y.^2.*x./r.^4 + x./r.^2);

Exx_x = dq_dx;
Eyy_x = dw_dy;
Exy_x = 0.5*(dw_dx + dq_dy);

Exx_y = dw_dx;
Eyy_y = dp_dy;
Exy_y = 0.5*(dp_dx + dw_dy);

%rotation rate
omega_x = .5*(dq_dy - dw_dx);
omega_y = .5*(dw_dy - dp_dx);


