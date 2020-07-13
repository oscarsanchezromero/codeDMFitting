function [aux]=rhoBUR(s)
% Numerical aproximation of the function rho_0(s) defined in Proposition 1.1 considering the Burkert model
%
% s is a row vector containing scale factors 
%
% rhoBUR(s) returns a row vector whose ith component is rho_0(s(i))
%
% This function requires global variables defined in the script redMethRotCurveFitting.m
global radii;
global vrot;
global weights;
global vbary;
global CteDim;
global somenullvbary;
global totalnullvbary;
%
aux=0*s;
vHalos=vBUR(radii,s);
rhs= WeighProd(vHalos,vHalos,weights); % rhs Eq 22 and 23 up to multiplicative constant 1/(s^3 CteDim) .
rhoVbaryNull= CteDim*(s'.^3).*(WeighProd(vrot*ones(size(s)),vHalos,weights)./rhs).^2;  % Eq 24
% This rootfinding can not be vectorizez because of the scalar character of arguments of routine fzero
% Those cases with vbary null at some radii have to be treated particularly.
if (totalnullvbary == true) % Case: vbary totally null
  aux = rhoVbaryNull';
elseif(somenullvbary == true) % Case: vbary null at some radii
    rango=find(rhs); % In this case condition  22 always holds.
    for i=rango'
      rhoequation=inline("rhs(i) - WeighProd(vrot,(vHalos(:,i).^2)./sqrt(t*(vHalos(:,i).^2)/((s(i)^3)*CteDim)+(vbary.^2)),weights)","t");
      j=-3;
      while (rhoequation(10^(j)*rhoVbaryNull(i)) > 0) j-- endwhile
      aux(i) = fzero(rhoequation,[10^(j)*rhoVbaryNull(i),rhoVbaryNull(i)]);
    endfor  
else  % Case: non null vbary at any radii 
  lhs=WeighProd(vrot*ones(size(s)),((vHalos).^2)./(vbary*ones(size(s))),weights);
  rango=find(rhs <= lhs); % checking condition  22
  for i=rango'
    rhoequation=inline("rhs(i) - WeighProd(vrot,(vHalos(:,i).^2)./sqrt(t*(vHalos(:,i).^2)/((s(i)^3)*CteDim)+(vbary.^2)),weights)","t");
    aux(i) = fzero(rhoequation,[0,rhoVbaryNull(i)]);
  endfor
endif
endfunction
