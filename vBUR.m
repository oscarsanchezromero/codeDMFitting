function [aux]=vBUR(x,s)
% Burkert rotation velocity (up to constants, see Eq.58).
%
% x is a column vector containing radii. 
% s is a row vector containing scale factors
%
% vBUR(x,s) returns a matrix where the (i,j) element 
% is the rotation velocity due to Burkert's density profile 
% at radius x(i) with scale factor s(j) and central density equal to 1
% up to multiplicative constants.
% See A. Burkert, Astrophys. J. 447 (1995) L25 for more details.
%
aux=sqrt(pi*(log((1+(x*s).^2).*((1+x*s).^2))-2*atan(x*s))./(x*ones(size(s))));
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
