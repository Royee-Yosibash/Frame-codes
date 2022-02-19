function [totErrorMeas] = getErrorMeasurement(dif, metricStr)
%GETERRORMEASUREMENT This function receives a vector of dif measurements
% and a metric to use return a single error value for that metric.
%
%	Input
% ------------------------
% (1) dif           -   A vector of diffrences/errors.
% (2) metricStr     -   A string that specifies the error measurement 
%                       needed. The current types supported are: 
%                       (a) L1
%                       (2) L2
% 
%   Output
% ------------------------
% (1) totErrorMeas	-	The error measurement value.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

if strcmpi(metricStr, 'L1')
    totErrorMeas = sum(abs(dif));
elseif strcmpi(metricStr, 'L2') || strcmpi(metricStr, 'MMSE')
    totErrorMeas = sum(dif.^2);
else
    error('No metric fitting the inputed metricStr');
end

end

