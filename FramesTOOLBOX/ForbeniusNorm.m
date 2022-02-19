function [norm] = ForbeniusNorm(mat)
%FORBENIUSNORM Calculates the Forbenius norm of an inputed matrice. The
%Forbenius norm is defined as the square root of the sumation over all the 
%absolute value and squared elements in the matrice

% Input
% -----------------------------
% (1) mat - A matrix of arbitrary size to calculate the Forbenius norm for.
%
% Output:
% -----------------------------
% (1) norm - The Forbenius norm for the matrix 'mat'.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021


siz = size(mat);

norm = 0;
for i = 1:siz(1)
    for j = 1:siz(2)
        norm = norm + abs(mat(i,j))^2;
    end
end
norm = sqrt(norm);

end

