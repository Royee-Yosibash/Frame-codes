function [f] = marchenkoPasturPDF(x, beta)
%MARCHENKOPASTURPDF This function receives a sample space vector and the
% parameter beta to return the Marchenko Pastur probability density
% function. An example for a description to this PDF can be found here: 
% https://en.wikipedia.org/wiki/Marchenko%E2%80%93Pastur_distribution
%
%	Input
% ------------------------
% (1) x     -   The span of the sample space. Inputed as a vector.
% (2) beta	-   The sub-frame aspect ratio.
% 
%   Output
% ------------------------
% (1) f     -	The Marchenko Pastur probability density function. Output 
%               is the same size as x.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

LambdaPlus = (1+sqrt(beta))^2;   LambdaMinus = (1-sqrt(beta))^2;

sqrfnc = (x >= LambdaMinus) & (x <= LambdaPlus);
f = sqrfnc .* sqrt( (x-LambdaMinus).*(LambdaPlus-x) ) ./ (2*pi()*beta*x);
f(x==0) = 0;

f = normPDF(x,f);

end

