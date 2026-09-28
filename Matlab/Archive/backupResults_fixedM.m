close all; clear all; clc
addpath(genpath('GUI'))

Code = [];
% castType = 'uint32';
%% Define Frame size
%   The system has 3 pramaters:
%   (1) n - number of workers
%   (2) m - number of workers if there were no stragglers (number of computations with no stragglers..)
%       (2.1)gamma = m/n
%   (3) k - number of computations returned
%       (3.1) beta = k/m

%   n>k>m

mComputations = 100;
P = 0.5; % = k/n = k/(m/gamma) =(k/m) * gamma = beta*gamma -> beta = P/gamma
gamma_vect = 0.3:0.02:P-0.02;

maxValue = 40;

%% Define Code and create initial non-normalized frame


codeTypes = {'Wishart'; 'LPF'; 'RS'};  
% codeTypes = {'BPF'; 'LPF' };
% codeTypes = {'LPF'};
% codeTypes = {'Wishart'; 'LPF'; 'RS'} ;
normDim = 'Coloumn';
numTests = 1000;

% check that if using a random set of the DFT we get manova
% check that the spectrum of RS is LPF?

%% Run iterations on number of returned nodes
plotResult = false;

%% Create the theoretical distributions's
t=0.001:0.001:10;
NoiseAmp_Vect = zeros(size(gamma_vect));
for jGamma = 1:numel(gamma_vect)
    gamma = gamma_vect(jGamma); nWorkers = round(mComputations/gamma);
    beta = P/gamma;
    [f] = manovaPDF(t, 1/beta, gamma);
    currEigen = getEigenValuesFromDist(t,f,mComputations);
    NoiseAmp_Vect(jGamma) = sum(1./currEigen)*sum(currEigen)/(mComputations^2);
end
NoiseAmp_dB = 10*log10(NoiseAmp_Vect);

hFigCompare = figure(); ax = axes(hFigCompare); hold(ax, 'on');
titleStr = {'Noise Amplification Vs. \gamma^{-1} for each code type'; ['Fixed k/n ratio of ', num2str(P) ]};
title(ax, titleStr);
xlabel(ax,'\gamma^{-1}');   ylabel(ax,'Noise Amplification [dB]');
plot(ax, (gamma_vect.^-1),NoiseAmp_dB, 'DisplayName', 'MANOVA\ETF');

%% Create Noise Amplification plots
for iCodeTypes =1:numel(codeTypes)    
    codeType = codeTypes{iCodeTypes};
    NoiseAmp_Vect = zeros(size(gamma_vect));
%     NoiseAmp_Vect = cast(NoiseAmp_Vect,castType);
    for jGamma = 1:numel(gamma_vect)
        gamma = gamma_vect(jGamma); nWorkers = round(mComputations/gamma);
        if mod(nWorkers-mComputations,2)
            nWorkers = nWorkers+1;
        end
        
        workers2Del= round((1-P) * nWorkers);
        NoiseAmp_Vect(jGamma) = collectFrameData(workers2Del, codeType, normDim, nWorkers, mComputations, gamma, numTests, plotResult);
    end
    if any(NoiseAmp_Vect<0)
        disp('Max gain value obtained, choosing maxvalue');
        disp([codeType, ', gamma= ' ,num2str(jGamma)]);
        NoiseAmp_Vect(NoiseAmp_Vect < 0) = 1.7976931348623158e+308-1; 
    end
    
    NoiseAmp_dB = 10*log10(abs(NoiseAmp_Vect));
    NoiseAmp_dB(NoiseAmp_dB>=maxValue) = maxValue;
    plot(ax, (gamma_vect.^-1),NoiseAmp_dB, '-o' ,'DisplayName', codeType);
end

grid(ax, 'on'); legend(ax);
% set(ax, 'XDir','reverse');
