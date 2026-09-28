%%

addpath(genpath('GUI'))

t = linspace(-10, 10, 10000);
beta = 0.5;
gamma = 0.5;
[fMP] = marchenkoPasturPDF(t, beta);
[fmanova] = manovaPDF(t, beta, gamma);

figure();
plot(t, fMP); hold on;
plot(t, fmanova);

%%
close all; clear all; clc
figure();   hold on;

numTests = 100;

%   The system has 3 pramaters:
%   (1) n - number of workers
%   (2) m - number of workers if there were no stragglers (number of computations with no stragglers..)
%   (3) k - number of computations returned
%   n>k>m

nWorkers = 160;
WorkerRedundancyNeeded  = 85 ;
nSets = nWorkers-WorkerRedundancyNeeded;

% nReturned = nWorkers - WorkerRedundancyNeeded + randi(WorkerRedundancyNeeded,1);
nReturned = nWorkers - WorkerRedundancyNeeded;

% Compare to DFT with some rows missing
DFTmat = (1/sqrt(nWorkers)) * dftmtx(nWorkers); % normalized DFT matrix

%% Randperm (maybe dif-set)
for i =1:numTests
    idx = sort(randperm(nWorkers,nReturned));
    DFTmatRedu = DFTmat(:,idx);
%     psuedoDFTInverseFilter = inv(DFTmatRedu' * DFTmatRedu) * DFTmatRedu';

    H = DFTmatRedu;
    A = transpose(H)*H;
    eigenValues = eig(A); % only for square...
    eigenValues_sqrt = sqrt(abs(eigenValues));
    unqEig = unique(eigenValues_sqrt);
    % A = H*transpose(H);
    % eigenValues = eig(A); % only for square...
    % eigenValues_sqrt = sqrt(abs(eigenValues))

    % s = svd(H);
    % [U,S,V] = svd(H);

    nBins = 2*numel(eigenValues_sqrt);
    maxX = max(abs(eigenValues_sqrt));
    x = linspace(-maxX, maxX, nBins);
    n = hist(eigenValues_sqrt, x);

    plot(x, n, 'b*'); 
end
%% LPF
% for i =1:numTests
    DFTmatLPF = DFTmat(:,1:numel(idx));
    % psuedoDFTInverseFilter = inv(DFTmatLPF' * DFTmatLPF) * DFTmatLPF';

    H = DFTmatLPF;
    A = transpose(H)*H;
    eigenValues = eig(A); % only for square...
    eigenValues_sqrt = sqrt(abs(eigenValues));
    unqEig = unique(eigenValues_sqrt);
    % A = H*transpose(H);
    % eigenValues = eig(A); % only for square...
    % eigenValues_sqrt = sqrt(abs(eigenValues))

    % s = svd(H);
    % [U,S,V] = svd(H);

    nBins = 2*numel(eigenValues_sqrt);
    maxX = max(abs(eigenValues_sqrt));
    x = linspace(-maxX, maxX, nBins);
    n = hist(eigenValues_sqrt, x);

    plot(x, n, 'r*'); 
% end

xlabel('eigenvalue'); ylabel('Count');

close all;
%%

figure(); hold on;
for i =1:(nWorkers-5)
    nReturned = nWorkers - i;
    idx = sort(randperm(nWorkers,nReturned));
    DFTmatLPF = DFTmat(:,1:numel(idx));
    H = DFTmatLPF;
    A = transpose(H)*H;
    eigenValues = eig(A); % only for square...
    eigenValues_sqrt = sqrt(abs(eigenValues));
    unqEig = unique(eigenValues_sqrt);
    nBins = 2*numel(eigenValues_sqrt);
    maxX = max(abs(eigenValues_sqrt));
    x = linspace(-maxX, maxX, nBins);
    n = hist(eigenValues_sqrt, x);
    plot(x, n, 'r*'); 
end
