function [f] = manovaPDF(x, beta, gamma)
%MANOVAPDF This function receives a sample space vector and the
% parameters beta and gamma to return the MAMNOVA density
% function. An example for a description to this PDF can be found in Royee 
% Yosibash's master thesis.
%
%
%	Input
% ------------------------
% (1) x     -   The span of the sample space. Inputed as a vector.
% (2) beta	-   The sub-frame aspect ratio.
% (3) gamma -   The frame aspect ratio.
% 
%   Output
% ------------------------
% (1) f     -	The MANOVA probability density function. Output is the same
%               size as x.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

if numel(x) == 1
    error('no span on x');
end

% Here beta is defined as m/k

leftTerm = sqrt((1-gamma)/beta);
rightTerm = sqrt(1-gamma/beta);

rPlus = (leftTerm + rightTerm)^2;   rMinus = (leftTerm - rightTerm)^2;

sqrfnc = (x >= rMinus) & (x <= rPlus);
aboveFrac = beta * sqrt((x-rMinus).*(rPlus-x));
belowFrac = (2*pi()*x).*(1-gamma*x);
f =  aboveFrac./belowFrac;
f = real(f) .* sqrfnc;
% f(t == (1/gamma)) = f(t == (1/gamma)) + max([0,(beta*gamma + gamma -1)])/min([beta*gamma , gamma]);

f = normPDF(x,f);

denom = gamma * ( 1 + 1/beta ) - 1;
if denom > 0
    deltaMult = denom/ min([gamma, gamma / beta]);
    dif = abs( x - (1/gamma) );
    [~,IndMin] = min(dif);
    
    f = f * (1-deltaMult); % Assuming it was normalized properly beforehand
    f(IndMin) = f(IndMin) + deltaMult; 
end



    
end

