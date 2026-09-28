function [matrices2decode, Acoded_partioned, Bcoded_partioned, EncodingMatA, EncodingMatB, mForCodeA, mForCodeB, isInPairs] = EncodingScheme(A, B, codeType, nNodes, mDataSets)
%%

[n1, N] = size(A);
[N, n2] = size(B);


isInPairs = false;
switch codeType
    case 'OrthoMatDot'
        mForCodeA = (mDataSets+1)/2;
        mForCodeB = mForCodeA;
    case 'Circulant Permutation'
        mForCodeA = mDataSets;
        mForCodeB = 1;
        
        isInPairs = true;
    otherwise
        mForCodeA = mDataSets;
        if size(B,2) == 1
            mForCodeB = 1;
        else      
            mForCodeB = mForCodeA;
        end
end

% Fix the matrices if they do not distribute equaly in M sub-matrices
[A_new, B_new] = FixMatricesDimensions(A, B, codeType, mForCodeA, mForCodeB);

%% 
CodeSubType = [];

isNormalize = false;
normDim = 'Column';
% normDim = 'None';
switch codeType
    case 'Circulant Permutation'
        [Code, isDeterministic, numPossibilities] = getCode(2*nNodes, 2*mForCodeA, codeType, CodeSubType);
    otherwise
        [Code, isDeterministic, numPossibilities] = getCode(nNodes, mForCodeA, codeType, CodeSubType);
end

EncodingMatA = transpose(Code);
if isNormalize
    switch normDim
        case 'Row'
            for iTemp = 1:size(EncodingMatA,2) % works on the inverse so I switched it
                EncodingMatA(:,iTemp) = EncodingMatA(:,iTemp)./norm(EncodingMatA(:,iTemp));
            end
            
        case 'Column'
            for iTemp = 1:size(EncodingMatA,1)
                EncodingMatA(iTemp,:) = EncodingMatA(iTemp,:)./norm(EncodingMatA(iTemp,:));
            end
        case 'None'
    end
end

EncodingMatB = EncodingMatA;

%% Encode
Acoded_partioned = [];
Bcoded_partioned = [];
matrices2decode = [];
switch codeType
    case 'OrthoMatDot'
        [Acoded_partioned] = encodeMatrix(A_new, EncodingMatA, 'col');
        [Bcoded_partioned] = encodeMatrix(B_new, EncodingMatB, 'row');
	case 'Circulant Permutation'
%         Acoded_partioned1 = encodeMatrix(A, EncodingMatA(:, 1:2:end), 'col');
%         Acoded_partioned2 = encodeMatrix(A, EncodingMatA(:, 2:2:end), 'col');
%         Acoded_partioned  = cell(size(Acoded_partioned1));
%         for i=1:numel(Acoded_partioned)
%             Acoded_partioned{i} = Acoded_partioned1{i}+Acoded_partioned2{i};
%         end
        Acoded_partioned = encodeMatrix(A_new, EncodingMatA, 'row');
        matrices2decode = cell(size(Acoded_partioned));
        for iTemp=1:numel(Acoded_partioned)
            matrices2decode{iTemp} = Acoded_partioned{iTemp} * B_new;
        end
%         clear Acoded_partioned
    otherwise
        [Acoded_partioned] = encodeMatrix(A_new, EncodingMatA, 'row');
%         [Bcoded_partioned] = encodeMatrix(B, EncodingMatB, 'col');
        matrices2decode = cell(size(Acoded_partioned));
        for iTemp=1:numel(Acoded_partioned)
            if size(B_new,2) == 1
                matrices2decode{iTemp} = Acoded_partioned{iTemp} * B_new(:,1);
            else
                error();
            end
        end
%         clear Acoded_partioned Bcoded_partioned
end


end

