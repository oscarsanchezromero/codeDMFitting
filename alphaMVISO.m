function [eval]=alphaMVISO(s)
% Evaluation of the function g(rho_0(s),s), where g is defined in equation 26
% considering the Isothermal model. By expression 19 this function provides
% varphi(s) adding appropriate terms. 
%
% s is a row vector containing scale factors 
%
% alphaMVISO(s) returns a column vector whose ith component is g(rho_0(s(i)),s(i))
%
% This function requires global variables defined in the script redMethRotCurveFitting.m
global radii;
global vrot;
global weights;
global vbary;
global CteDim;
%
rhoaux=(rhoISO(s));
vaux=vISO(radii,s);
eval=rhoaux'.*WeighProd(vaux,vaux,weights)./(CteDim*s'.^3);
eval=eval-2*(WeighProd(vrot*ones(size(s)),sqrt(((vaux).^2).*(ones(length(radii),1)*(rhoaux./(CteDim*s.^3)))+(vbary.^2)*ones(size(s))),weights));
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

