function [eigenValues, conditionNumber] = getGramMatrixEigenvalues(H)
%GETGRAMMATRIXEIGENVALUES Calculates the gram matrix eigenvalues of some
%matrix H. If desired the condition number is also an available output.
%
% Input
% -------------------------
% (1) H                 -   Some matrix of arbitrary size of size [m,n]. 
%
% Output
% -------------------------
% (1) eigenValues       -	A vector of eignvaluess of length min([n,m])
% (2) conditionNumber	-	A double that represents the condition number  
%                           of the gram matrix. 
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

if size(H,1) > size(H,2)
    A = ctranspose(H)*H;
else
    A = H*ctranspose(H);
end

eigenValues = (eig(A));
conditionNumber = cond(A);
end