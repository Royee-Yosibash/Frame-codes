function [Eigen] = getEigenValuesFromDist(X,PDF,numEigenValues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CDF = cumsum(PDF);
R = rand(numEigenValues,1);
EigenIdx = zeros(1,numel(R));
for i=1:numel(R)
    EigenIdx(i)=find(R(i)<CDF,1,'first');
end
Eigen = X(EigenIdx);

end

