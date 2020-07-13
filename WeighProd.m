function [prod]=WeighProd(x,y,sigmas)
% Weighted scalar product function.
% x and y are matrices of order N times a and N times b, whose columns 
% are the factors for the weighted product.
% sigmas is a vector of lenght N containing the weights. 
%
% WeighProd(x,y,sigmas) returns a vector where the ith element 
% is the weighted scalar product of the ith column of x times the 
% ith column of y. i ranges between 1 and the minimum{a,b}
% 
prod= diag(x'*diag(sigmas)*y);
end 
