close all; clear all; clc
Code = [];
%% Define Frame size
%   The system has 3 pramaters:
%   (1) n - number of workers
%   (2) m - number of workers if there were no stragglers (number of computations with no stragglers..)
%       (2.1)gamma = m/n
%   (3) k - number of computations returned
%       (3.1) beta = k/m

%   n>k>m
nWorkers = 160;
WorkerRedundancyNeeded  = 80;
mComputations = nWorkers - WorkerRedundancyNeeded;
gamma = mComputations/nWorkers;

workers2Del=0:2:(WorkerRedundancyNeeded-20);
%% Define Code and create initial non-normalized frame
% codeType = 'RS';    normDim = 'Coloumn'; numTests = 100;
% codeType = 'Wishart';    normDim = 'None';  numTests = 2000;
% codeType = 'Wishart';    normDim = 'Coloumn';  numTests = 2000;
% codeType = 'DFT';    normDim = 'Coloumn'; numTests = 100;

% codeTypes = {'Wishart'; 'RS'; 'LPF' };  
codeTypes = {'BPF'; 'LPF' };
% codeTypes = {'Wishart'; 'LPF'; 'RS'} ;
normDim = 'Coloumn';
numTests = 1000;

% check that the spectrum of the DFT and RS are LowPass
% check that if using a random set of the DFT we get manova
% check that the spectrum of RS is LPF?

%% Run iterations on number of returned nodes
plotResult = false;
hFigCompare = figure(); ax = axes(hFigCompare); hold(ax, 'on');
for i =1:numel(codeTypes)
    codeType = codeTypes{i};
    [NoiseAmp_Vect] = collectFrameData(workers2Del, codeType, normDim, nWorkers, mComputations, gamma, numTests, plotResult);
    NoiseAmp_dB = 10*log10(NoiseAmp_Vect);
    plot(ax, 1-(workers2Del/nWorkers),NoiseAmp_dB, 'DisplayName', codeType);
end

grid(ax, 'on'); legend(ax);
xlabel(ax,'Workers returns / number of workers')
ylabel(ax,'Noise Amplification [dB]');
set(ax, 'XDir','reverse');
