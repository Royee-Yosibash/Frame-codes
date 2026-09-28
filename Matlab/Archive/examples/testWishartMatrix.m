close all; clear all; clc;
addpath(genpath('./FramesTOOLBOX'))

n = 400;
m = 200;
k = n;
gamma = m/n;
beta = k/m;
var = 1;

numTests = 10000;
numEigen = m;

eigenMat = inf*ones(numTests, numEigen);
for i = 1:numTests
    H = 0 + sqrt(var).*randn(m,n);
    A = (1/n) * H * transpose(H);
    eigenValues = eig(A); % only for square...
%     eigenValues_sqrt = sqrt(abs(eigenValues));
%     eigenMat(i, :) = eigenValues_sqrt;
    eigenMat(i, :) = eigenValues;
end
eigenVect = reshape(eigenMat, [1, numel(eigenMat)]);

nBins = 20*numel(eigenVect)/numTests;

titleStr = ['eignvalue probabilty for m=', num2str(m), ' n=', num2str(n)];
[hFig] = plotFramePDF(eigenVect, nBins, gamma, gamma, titleStr);
% [hFig] = plotFramePDF(eigenVect, nBins, beta, gamma, titleStr);
