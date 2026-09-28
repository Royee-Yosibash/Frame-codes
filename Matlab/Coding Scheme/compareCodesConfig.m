function config = compareCodesConfig()
%COMPARECODESCONFIG Default configuration for compareCodes.
%   The values here match the original defaults in compareCodes.m.

config.numBits = 8;
config.mDataSets = [29];
config.nNodes = [31];
config.SNR = [80];
config.numTrials = 1*10^4;
config.N = 2800;
config.n2 = 1;
config.codeTypes = {'Non consecutive powers'; 'Circulant Permutation'};

config.measuresRealName = {'MSE'; 'Frobenius norm'; ...
    'mean(conditionNumber)'; 'min(conditionNumber)'; 'max(conditionNumber)'};
config.measuresLatex = {'MSE'; '$\frac{\|\widehat{\mathbf{A}\mathbf{x}} - \mathbf{A}\mathbf{x}\|_F}{\|\mathbf{A}\mathbf{x}\|_F}$'; ...
    '$mean(\kappa_{\mathbf{E}_{dec}})$'; '$min(\kappa_{\mathbf{E}_{dec}})$'; '$max(\kappa_{\mathbf{E}_{dec}})$'};

end
