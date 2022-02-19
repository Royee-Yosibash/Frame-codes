function [newf] = normPDF(x,f)
%NORMPDF This function receives a probability density function described by
% the span of the sample space and the likelyhood value and returns a new
% (or equal) likelihood value so that the summation of the likelyhood
% values sum to 1.
%	Input
% ------------------------
% (1) x     -   The span of the sample space. Inputed as a vector.
% (2) f     -   The likelyhood value of the sample space. Inputed as a 
%               in the same size as x.
% 
%   Output
% ------------------------
% (1) newf	-	The new value of the probability density function at the
%               same sample point after normalizing. Output is the same  
%               size as x.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

dx = abs(x(2)-x(1));
sumUnderCurve = sum(f*dx);
newf = f/sumUnderCurve;

end

