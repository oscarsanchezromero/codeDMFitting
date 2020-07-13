function [aux]=vISO(x,s)
% Isothermal rotation velocity (up to constants, see Eq.58).
%
% x is a column vector containing radii. 
% s is a row vector containing scale factors
%
% vISO(x,s) returns a matrix where the (i,j) element 
% is the rotation velocity due to Isothemal density profile 
% at radius x(i) with scale factor s(j) and central density equal to 1
% up to multiplicative constants.
% See J. Gunn, J.R. Gott, Astrophys. J. 176 (1972) for more details.
% 
aux=sqrt(4*pi*(x*s-atan(x*s))./(x*ones(size(s))));
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
