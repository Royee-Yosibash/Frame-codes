clc; clear all;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup();
cd(fullfile(projectRoot, 'Coding Scheme'));
%% Create results table
config = compareCodesConfig();
measuresRealName = config.measuresRealName;
measuresLatex = config.measuresLatex;

%% Machine features/properties
numBits = config.numBits; % 8 bits
Imax = 2^numBits;
mDataSets = config.mDataSets;
nNodes = config.nNodes;

% mDataSets = 80;
% nNodes = 100;


% beta = mDataSets/(mDataSets+RedundancyLeft);
% betaInv = 1/beta;
SNR = config.SNR; % dB
numTrials = config.numTrials;
%% Matrices to be calculated
% A = randi(Imax, n1,N) - 1;
% B = randi(Imax, N, n2) - 1;
N = config.N; % A and B's common Matrix dimensions
n2 = config.n2; % B's number of columns

% mu = 0; sigma = 1;
% A = normrnd(mu, sigma, n1, N);
% B = normrnd(mu, sigma,  N, n2);

%% Simulation Inputs
% codeType = 'OrthoMatDot'
% codeType = 'Omitted Vandermonde'
codeTypes = config.codeTypes;

%% Check inputs
if numel(nNodes) ~= numel(mDataSets)
    error('number of nodes and datasets mismatch');
elseif ~all((nNodes-mDataSets)>=0)
    error('number of nodes must be greater or equal than the number of datasets');
end
%% Run coding schemes

TotResults = {};
for kSNR = 1:numel(SNR)
    for iNM = 1:numel(nNodes)
        disp(['Number of workers ', num2str(nNodes(iNM)), '. Number of original workers ', num2str(mDataSets(iNM))])
        newResultsStruct = SingleCodingSchemeResults(nNodes(iNM), mDataSets(iNM), SNR(kSNR), measuresRealName, measuresLatex);
        
        n1 = 22*mDataSets(iNM)*30; % A's number of rows
        A = randi(101, n1, N) - 51; B = randi(101,  N, n2) - 51;
        C = A*B;    % Real Result
        NumStragglers = 0:1:(nNodes(iNM)-mDataSets(iNM));
        for iCode = 1:numel(codeTypes)
            [matrices2decode, Acoded_partioned, Bcoded_partioned, EncodingMatA, EncodingMatB, mForCodeA, mForCodeB, isInPairs] = ...
                EncodingScheme(A, B, codeTypes{iCode}, nNodes(iNM), mDataSets(iNM));
            
            [meanCondNumber,maxCondNumber, minCondNumber] = CalculateGMConditionNumberValues(EncodingMatA, numTrials, NumStragglers, isInPairs);
            
            for jStraggler=1:numel(NumStragglers)
                disp(['Number of nodes left: ', num2str(mDataSets(iNM)+NumStragglers(jStraggler))]);
                k = nNodes(iNM)-NumStragglers(jStraggler);
                [MSE_mean, ForbeniusNormError_mean] = ...
                    DecodingScheme(C, codeTypes{iCode}, matrices2decode, EncodingMatA, EncodingMatB, nNodes(iNM), mDataSets(iNM), mForCodeA, mForCodeB, NumStragglers(jStraggler), SNR(kSNR), numTrials);
                %                     Colname = [codeTypes{iCode}, ' k is ', num2str(k)];
                newResultsStruct.addTableCol(k, codeTypes{iCode}, {MSE_mean; ForbeniusNormError_mean; ...
                    meanCondNumber(jStraggler); minCondNumber(jStraggler); maxCondNumber(jStraggler)});
            end
        end
        newResultsStruct.writeToTXTfile();
        TotResults = [TotResults, newResultsStruct];
        
    end
end

plotAllResults;
