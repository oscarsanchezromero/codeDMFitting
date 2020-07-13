function [aux]=vNFW(x,s)
% Navarro-Frenk-White rotation velocity (up to constants, see Eq.58).
%
% x is a column vector containing radii. 
% s is a row vector containing scale factors
%
% vISO(x,s) returns a matrix where the (i,j) element 
% is the rotation velocity due to Isothemal density profile 
% at radius x(i) with scale factor s(j) and central density equal to 1
% up to multiplicative constants.
% See Navarro, J. F., Frenk, C. S., & White, S. D. M. 1996, ApJ, 462, 563 for more details.
% 
aux=sqrt(4*pi*(log(1+x*s)./(x*ones(size(s)))-(ones(size(x))*s)./(1+x*s)));
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
