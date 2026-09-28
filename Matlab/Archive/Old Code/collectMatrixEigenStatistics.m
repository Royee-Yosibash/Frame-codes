function [eigenVect,NoiseAmp] = collectMatrixEigenStatistics(codeType, normDim, n, m, k, numTests)

% numEigen = n - m;
numEigen = m;
eigenMat = inf * ones(numTests, numEigen);
NoiseAmp = zeros(numTests, 1);
H = [];
for j = 1:numTests
    if isempty(H) 
        [H, isDeterministic] = getCode(n, m, codeType, normDim);
    elseif ~isDeterministic
        [H, isDeterministic] = getCode(n, m, codeType, normDim);
    end
    
    idx = sort(randperm(n,k));
    if size(H,1) > size(H,2)
        H_idx = H(idx,:);
    else
        H_idx = H(:, idx);
    end
    
    [eigenValues] = getGramMatrixEigenvalues(H_idx);
    eigenMat(j, :) = eigenValues;
    NoiseAmp(j) = sum(1./eigenMat(j, :))*sum(eigenMat(j, :))/(m^2);
end

eigenVect = reshape(eigenMat, [1, numel(eigenMat)]);
end

