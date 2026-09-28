close all; clear all; clc;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup();
cd(projectRoot);

%% MANOVA dist

titleStr = 'The MANOVA probability density function';

gamma = 0.5;
beta = 0.8;     % here beta is m/k
% [lambda,PDF] = createPDFOfEigenvalues(eigenVect,nBins);
% plot(lambda, PDF, 'b*', 'DisplayName','Gram matrix \lambda distribution'); grid on;


% meanEig = sum(lambda.*PDF);
% plot(meanEig, 0, 'bx', 'markersize', 20, 'DisplayName','Mean of Gram matrix \lambda');
% plot(beta, 0, 'g+', 'markersize', 20, 'DisplayName','\beta');

lambdaVect = 0:0.0001:10;

[f] = manovaPDF(lambdaVect, beta, gamma); 
hFig = figure('units','normalized','outerposition',[0 0 1 1]);   hold on;
ylabel('Probability Density'); xlabel('\lambda');
plot(lambdaVect, f, 'k', 'DisplayName','MANOVA PDF');
% plot(meanF, 0, 'rx', 'markersize', 20, 'DisplayName','Mean of marchenko pastur dist');
ymaxPlot = find(f>1e-8);  ymaxPlot = ymaxPlot(end);
xlim([0, 1.3*lambdaVect(ymaxPlot)]);
title(titleStr);    
% legend();


%% MANOVA example with DFT

titleStr = 'The MANOVA probability density function';

n = 400;
m = n*gamma;
k = m/beta;
NumStragglers = n-k;

Code = getCode(n, m, 'Omitted Vandermonde', 'None');

for i=1:size(Code, 2)
    Code(:,i) = Code(:,i)/norm(Code(:,i));
end
numTests = 1000;
numEigen = m;

eigenMat = inf*ones(numTests, numEigen);
for i = 1:numTests
    
    idxALeft = randperm(n);
    idxALeft = sort(idxALeft(NumStragglers+1:end));
    H = Code(:,idxALeft);
    [eigenValues] = getGramMatrixEigenvalues(H);
%     eigenValues_sqrt = sqrt(abs(eigenValues));
%     eigenMat(i, :) = eigenValues_sqrt;
    eigenMat(i, :) = eigenValues;
end
eigenVect = reshape(eigenMat, [1, numel(eigenMat)]);

nBins = 20*numel(eigenVect)/numTests;

titleStr = ['eignvalue probabilty for m=', num2str(m), ' n=', num2str(n)];



denom = gamma * ( 1 + 1/beta ) - 1;
if denom>0
    [lambdaVect,PDF] = createPDFOfEigenvalues(eigenVect,nBins,gamma);
else
    [lambdaVect,PDF] = createPDFOfEigenvalues(eigenVect,nBins,[]);
end

[f] = manovaPDF(lambdaVect, beta, gamma); 
hFig = figure('units','normalized','outerposition',[0 0 1 1]);   hold on;
ylabel('Probability Density'); xlabel('\lambda');
plot(lambdaVect, f, 'r', 'DisplayName','MANOVA PDF');
% plot(meanF, 0, 'rx', 'markersize', 20, 'DisplayName','Mean of marchenko pastur dist');
plot(lambdaVect, PDF, 'k*', 'DisplayName','Gram matrix \lambda hist'); grid on;


%% MP dist
clear all;
titleStr = 'The Marchenko Pastur probability density function';

gamma = 0.5;
beta = 0.8;
% [lambda,PDF] = createPDFOfEigenvalues(eigenVect,nBins);
% plot(lambda, PDF, 'b*', 'DisplayName','Gram matrix \lambda distribution'); grid on;


% meanEig = sum(lambda.*PDF);
% plot(meanEig, 0, 'bx', 'markersize', 20, 'DisplayName','Mean of Gram matrix \lambda');
% plot(beta, 0, 'g+', 'markersize', 20, 'DisplayName','\beta');

lambdaVect = 0:0.0001:10;

[f] = marchenkoPasturPDF(beta, lambdaVect); % This is originaly beta, but due to gamma being defined as m/n this is the one

% meanF = sum(lambda.*f);
hFig = figure('units','normalized','outerposition',[0 0 1 1]);   hold on;
ylabel('Probability Density'); xlabel('\lambda');
plot(lambdaVect, f, 'k', 'DisplayName','Marchenko Pastur (MP) PDF');
% plot(meanF, 0, 'rx', 'markersize', 20, 'DisplayName','Mean of marchenko pastur dist');
ymaxPlot = find(f>1e-8);  ymaxPlot = ymaxPlot(end);
xlim([0, lambdaVect(ymaxPlot)]);
title(titleStr);    
% legend();


%% Show for Random wishart matrix
n = 400;
m = n*beta;
var = 1;

numTests = 10000;
numEigen = m;

eigenMat = inf*ones(numTests, numEigen);
for i = 1:numTests
    H = 0 + sqrt(var).*randn(m,n);
    H = sqrt(1/n) * H;
    [eigenValues] = getGramMatrixEigenvalues(H);
%     eigenValues_sqrt = sqrt(abs(eigenValues));
%     eigenMat(i, :) = eigenValues_sqrt;
    eigenMat(i, :) = eigenValues;
end
eigenVect = reshape(eigenMat, [1, numel(eigenMat)]);

nBins = 20*numel(eigenVect)/numTests;

titleStr = ['eignvalue probabilty for m=', num2str(m), ' n=', num2str(n)];


[lambdaVect,PDF] = createPDFOfEigenvalues(eigenVect,nBins, []);

[f] = marchenkoPasturPDF(beta, lambdaVect);
hFig = figure('units','normalized','outerposition',[0 0 1 1]);   hold on;
ylabel('Probability Density'); xlabel('\lambda');
plot(lambdaVect, f, 'r', 'DisplayName','Marchenko Pastur (MP) PDF', 'linewidth', 1);
plot(lambdaVect, PDF, 'k*', 'DisplayName','Gram matrix \lambda hist'); grid on;
