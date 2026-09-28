function [B] = colSwitch(A, i ,j)
%COLSWITCH Summary of this function goes here
%   Detailed explanation goes here

n = size(A,2);
I = eye(n);
Inew = I;

Inew(:,i) = I(:,j); 
Inew(:,j) = I(:,i); 

B = A * Inew;

end

