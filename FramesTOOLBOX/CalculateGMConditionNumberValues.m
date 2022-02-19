function [meanCondNumber,maxCondNumber, minCondNumber] = CalculateGMConditionNumberValues(inputMat, numTrials, num2Remove, isInPairs)
%CALCULATEGMCONDITIONNUMBERVALUES This function receives some matrice 
%(inputMat) and calculates the matrix's gram matrix's condition number as
% many times as required by the input and with regard to the number of
% straglers.
%
%	Input
% ------------------------
% (1) inputMat      -   The matrix to calculate the condition number for.
% (2) numTrials     -   The number of trails to go over and calculate the
%                       condition number for.
% (3) num2Remove	-   The number of rows that needs to be randomly
%                       deleted in each iteration.
% (4) isInPairs     -   A special flag that shoud be o be turned on if two
%                       consecutive calculations are done in the same node.
% 
%   Output
% ------------------------
% (1) meanCondNumber	-	The mean of the condition numbers calculated
%                           over all the iterations.
% (2) maxnCondNumber	-	The maximal condition number calculated over
%                           all the iterations.
% (3) minnCondNumber	-	The minimal condition number calculated over
%                           all the iterations.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021



nNodes = size(inputMat,1);

meanCondNumber = zeros(size(num2Remove));
maxCondNumber = meanCondNumber;
minCondNumber = meanCondNumber;
for jStraggler = 1:numel(num2Remove)
    
    conditionNumber = zeros(1,numTrials);
    NumStraggler = num2Remove(jStraggler);
    for iTrial=1:numTrials
        if ~isInPairs
        idxALeft = randperm(nNodes);
        idxALeft = sort(idxALeft(NumStraggler+1:end));
        else
            idxALeft = randperm(nNodes/2);
            idxALeft = sort(idxALeft(NumStraggler+1:end));
            idxALeft = sort([2*idxALeft-1, 2*idxALeft]);
        end
        LeftOverMat = inputMat(idxALeft,:);
        if size(LeftOverMat,1) == size(LeftOverMat,2)
            decoder = inv(LeftOverMat);
            eigenValues = abs(eig(decoder));
            conditionNumber(iTrial) = max(eigenValues)/min(eigenValues);
        else
            decoder = inv(ctranspose(LeftOverMat)* LeftOverMat) * ctranspose(LeftOverMat);
%             [~,e,~] = svd(decoder);
%             e = diag(e);
%             e(e==0) = [];
%             e = abs(e);
            [~, conditionNumber_Gram] = getGramMatrixEigenvalues(decoder);
            conditionNumber(iTrial) = sqrt(conditionNumber_Gram);
        end

        
    end

    meanCondNumber(jStraggler) = mean(conditionNumber);
    maxCondNumber(jStraggler) = max(conditionNumber);
    minCondNumber(jStraggler) = min(conditionNumber);
end

end

