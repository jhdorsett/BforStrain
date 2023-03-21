function [u,v,exx] = EdgeDisp_finite(x_obs, disloc)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Displacements Due to Two Edge Dislocations
% (fault has finite length down dip)
%
% Paul Segall 1996
% Input:
%		x_obs = coordinates of the observation points
%		disloc(1) = depth of updip end;
%		disloc(2) = horizontal position of updip end
%		disloc(3) = length downdip;
%		disloc(4) = dip (degrees);
%		disloc(5) = slip;
% modification by Kaj Johnson, May 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dip = disloc(4)*pi/180;
L   = disloc(3);
d1  = disloc(1);
d2  = disloc(1) + L*sin(dip);
s_v  = disloc(5)*sin(dip);
s_h  = disloc(5)*cos(dip);
one = ones(size(x_obs)) ;

x_obs=x_obs-disloc(2); %shift for horizontal position of updip end

zeta_1 = x_obs / d1;
zeta_2 = (x_obs - L*cos(dip)*one) /d2;
denom_1  = one + zeta_1.^2;
denom_2  = one + zeta_2.^2;

%vertical displacement
v = ( s_v*atan(zeta_1) + (s_h*one + s_v*zeta_1)./denom_1 ...
 	-s_v*atan(zeta_2) - (s_h*one + s_v*zeta_2)./denom_2 )/pi;

%horizontal displacement
u = -( s_h*atan(zeta_1) + (s_v*one - s_h*zeta_1)./denom_1 ...
 	-s_h*atan(zeta_2) - (s_v*one - s_h*zeta_2)./denom_2 )/pi;

% horizontal strain -NOTE normalized by d2
exx = (d2/d1)*(s_v*zeta_1 - s_h*zeta_1.^2)./ denom_1.^2 ...
	 -(s_v*zeta_2 - s_h*zeta_2.^2)./ denom_2.^2;
exx = 2*exx/pi;
