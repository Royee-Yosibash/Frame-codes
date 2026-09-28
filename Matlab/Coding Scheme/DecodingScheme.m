function [MSE_mean, ForbeniusNormError_mean] = DecodingScheme(C, codeType, matrices2decode, EncodingMatA, EncodingMatB, nNodes, mDataSets, mForCodeA, mForCodeB, NumStragglers, SNR, numTrials)
%DECODINGSCHEME Summary of this function goes here
%   Detailed explanation goes here

MSE_mean = 0;   ForbeniusNormError_mean = MSE_mean;
[n1,n2] = size(C);
dynamicRangeFactor = 10^29; % this is just so calculations dont go over the bit limit

%% Run decoding trials
noisePower = 10^(-SNR/10);
MSE = zeros(1,numTrials);
meanError = zeros(1,numTrials);
ForbeniusNormError = zeros(1,numTrials);

for iTrial = 1:numTrials
    if mod(iTrial,1000) == 0
        disp(['Trial number ', num2str(iTrial), '#']);
    end
    
    %% Add Stragglers
    idxALeft = randperm(nNodes);
    idxALeft = sort(idxALeft(NumStragglers+1:end));
    
    %% Decode the soloution
    C_reconstructed = [];
    switch codeType
        case 'OrthoMatDot'
            if numel(idxALeft) ~= mDataSets
                error('This decoder can not use more nodes then necessary');
            end
            
            matrices2decodeLeft = matrices2decode(idxALeft);
            for iTemp=1:numel(matrices2decodeLeft)
                matrices2decodeLeft{iTemp} = matrices2decodeLeft{iTemp} +randn(size(matrices2decodeLeft{iTemp})).* sqrt(noisePower);
            end
            
            C_forReconstruction = matrices2decodeLeft;
            
            [G_P, ~, ~] = getCode(nNodes, (2*mForCodeA-1),  codeType, CodeSubType);
            G_R = G_P(:,idxALeft);
            
            % Fusion node procedure from pseudo-code
            Ginv = inv(G_R);
            C_reconstructed = zeros(n1,n2);
            
            [G_second, ~, ~] = getCode(mForCodeA, (2*mForCodeA-1), codeType, CodeSubType);
            
            for iTemp = 1:n1
                for jTemp = 1:n2
                    
                    coeff = zeros(1,(2*mForCodeA-1));   vectAnswers = coeff;
                    for jj = 1:numel(vectAnswers)
                        temp = C_forReconstruction{jj};
                        vectAnswers(jj) = temp(iTemp,jTemp);
                    end
                    
                    % start line 1
                    coeff = vectAnswers * Ginv;
                    % end line 1
                    
                    % start line 2
                    v = coeff*G_second; % is this the same vector as vectAnswers?
                    % end line 2
                    
                    % start line 3
                    C_reconstructed(iTemp,jTemp) = (2/mForCodeA)*sum(v);
                    % end line 3
                    C_reconstructed(iTemp,jTemp) = C_reconstructed(iTemp,jTemp); % for m=4 it was 16????
                    
                end
            end
            % End psuedocode
            
            ratio = C_reconstructed./C;
            meanRatio = mean(mean(ratio))
            C_reconstructed = C_reconstructed./meanRatio;
            
        case 'Circulant Permutation'
            idxALeft = sort([2*idxALeft-1, 2*idxALeft]);
            matrices2decodeLeft = matrices2decode(idxALeft);
            for iTemp=1:numel(matrices2decodeLeft)
                matrices2decodeLeft{iTemp} = matrices2decodeLeft{iTemp} +randn(size(matrices2decodeLeft{iTemp})).* sqrt(noisePower);
            end
            
            % decode using pseudocode
            LeftOverEncodingMat = EncodingMatA(idxALeft,:); %Grot
            
            if size(LeftOverEncodingMat,1) == size(LeftOverEncodingMat,2)
                decoder = inv(LeftOverEncodingMat);
            else
                decoder = inv(ctranspose(LeftOverEncodingMat)* LeftOverEncodingMat) * ctranspose(LeftOverEncodingMat);
            end
            
            
            C_decoded = cell(1,2*mDataSets);
            for iTemp = 1:numel(C_decoded)
                decodingRow = decoder(iTemp,:);
                C_decoded{iTemp}= 0;
                for jTemp = 1:length(decodingRow)
                    C_decoded{iTemp} = C_decoded{iTemp} + decodingRow(jTemp) * matrices2decodeLeft{jTemp}; % change
                end
            end
            
            C_reconstructed = [];
            for iTemp = 1:length(C_decoded)
                C_reconstructed = [C_reconstructed ;  C_decoded{iTemp}];
            end
            
            
        case {'Omitted Vandermonde', 'Non consecutive powers'}
            matrices2decodeLeft = matrices2decode(idxALeft);
            for iTemp=1:numel(matrices2decodeLeft)
                matrices2decodeLeft{iTemp} = matrices2decodeLeft{iTemp} +randn(size(matrices2decodeLeft{iTemp})).* sqrt(noisePower);
            end
            
            % decode using pseudocode
            LeftOverEncodingMat = EncodingMatA(idxALeft,:); %Grot
            
            if size(LeftOverEncodingMat,1) == size(LeftOverEncodingMat,2)
                decoder = inv(LeftOverEncodingMat);
            else
                decoder = inv(ctranspose(LeftOverEncodingMat)* LeftOverEncodingMat) * ctranspose(LeftOverEncodingMat);
            end
            
            C_decoded = cell(1,mDataSets);
            for iTemp = 1:numel(C_decoded)
                decodingRow = decoder(iTemp,:);
                C_decoded{iTemp}= 0;
                for jTemp = 1:length(decodingRow)
                    C_decoded{iTemp} = C_decoded{iTemp} + decodingRow(jTemp) * matrices2decodeLeft{jTemp}; % change
                end
            end
            
            C_reconstructed = [];
            for iTemp = 1:length(C_decoded)
                C_reconstructed = [C_reconstructed ;  C_decoded{iTemp}];
            end
            
            
    end
    
    C_reconstructed = C_reconstructed(1:size(C,1), 1:size(C,2));
    difC = C_reconstructed - C;
    difPrecC = 100*abs(difC)./C;
    
    meanError(iTrial) = mean(mean(difPrecC(1:n1,1:n2)))*10^29;
    
    ForbeniusNormError(iTrial) = ForbeniusNorm(difC)/ForbeniusNorm(C)*dynamicRangeFactor;
    MSE(iTrial) = sum(abs((difC).^2))/sum(abs(C_reconstructed.^2))*dynamicRangeFactor;
end

MSE_mean = mean(MSE)/dynamicRangeFactor;
ForbeniusNormError_mean= mean(ForbeniusNormError)/dynamicRangeFactor; %forbenius norm
% mean_meanError= mean(meanError)*10^-29

end

