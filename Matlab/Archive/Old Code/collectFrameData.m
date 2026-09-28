function [NoiseAmp_Vect] = collectFrameData(workers2Del, codeType, normDim, nWorkers, mComputations, gamma, numTests, plotResult)
%COLLECTFRAMEDATA Summary of this function goes here
%   Detailed explanation goes here
NoiseAmp_Vect = zeros(size(workers2Del));
for i =1:numel(workers2Del)
    k = nWorkers - workers2Del(i); % nReturned
    beta = k/mComputations;
%     if beta > 1
%         beta = beta^-1;
%     end
    
    [eigenVect,NoiseAmp] = collectMatrixEigenStatistics(codeType, normDim, nWorkers, mComputations, k, numTests);    
    NoiseAmp_Vect(i) = mean(NoiseAmp);
    
    if plotResult
        titleStr = {['eignvalue probabilty for \beta = ' , num2str(beta) , ' and \gamma = ', num2str(gamma)], ...
            [num2str(k), ' workers returned (=k) out of ', num2str(nWorkers), ' initiated (=n). ', num2str(mComputations), ' Computations are needed (=m)'], ...
            ['Noise Amplification = ', num2str(NoiseAmp_Vect(i))]};
        nBins = 20*numel(eigenVect)/numTests;
        [hFig] = plotFramePDF(eigenVect, nBins, beta, gamma, titleStr);
    end

    disp(['The noise amplication is: ', num2str(NoiseAmp_Vect(i))]);    
end

end

