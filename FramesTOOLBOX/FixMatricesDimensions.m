function [Anew, Bnew] = FixMatricesDimensions(A, B, codeType, mForCodeA, mForCodeB)
%FIXMATRIXDIMENSIONS This function receives a two matrices, A and B, that
% should be multipled by dot product as A*B and fixed their dimension to 
% fit the encoding process. 
%
% Explanation: If the matrices A and B are to be encoded with some code the 
% first step is to partition them to some mForCodeA and mForCodeB parts 
% respectivly. If the columns/rows of these matrices do not fit the
% partition in a way that creates equal partitioned matrices that are all
% equaly sized then columns/rows of zeros are inserted to create an equal
% partition of the columns/rows. These columns/rows in the final dot  
% product A*B will produce an expected number of rows that the decoder can
% omit.
%
%	Input
% ------------------------
% (1) A             -   The matrix A from the product A*B.
% (2) B             -   The matrix B from the product A*B.
% (3) codeType      -   The type of code that is used for
%                       encoding/decoding.
% (4) mForCodeA     -   The number of partitions needed for encoding A.
% (5) mForCodeB     -   The number of partitions needed for encoding B.
%
%   Output
% ------------------------
% (1) Anew          -	The matrix A from the product A*B after fixing
%                       dimensions.
% (2) Bnew          -	The matrix B from the product A*B after fixing
%                       dimensions.
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

sizA = size(A);
sizB = size(B);

Anew = A;   Bnew = B;

switch codeType
    case 'OrthoMatDot'
        if mod(sizA(2),mForCodeA) ~= 0
            Anew = [A, zeros(sizA(1), mForCodeA - mod(sizA(2),mForCodeA))];
            Bnew = [B; zeros(mForCodeB - mod(sizB(1),mForCodeB), sizB(2))];
        end
        
    case 'Circulant Permutation'
%         if mod(sizA(2),mForCodeA) ~= 0
%             Anew = [A, zeros(sizA(1), mForCodeA - mod(sizA(2),mForCodeA))];
%         end
%         
%         Bnew = [B; zeros(mForCodeB - mod(sizB(1),mForCodeB), sizB(2))];
        
        if mod(sizA(1),mForCodeA) ~= 0
            Anew = [A; zeros( mForCodeA - mod(sizA(1),mForCodeA), sizA(2))];
        end

        if mod(sizB(2),mForCodeB) ~= 0
            Bnew = [B, zeros(sizB(1), mForCodeB - mod(sizB(2),mForCodeB))];
        end

    otherwise
        if mod(sizA(1),mForCodeA) ~= 0
            Anew = [A; zeros( mForCodeA - mod(sizA(1),mForCodeA), sizA(2))];
        end

        if mod(sizB(2),mForCodeB) ~= 0
            Bnew = [B, zeros(sizB(1), mForCodeB - mod(sizB(2),mForCodeB))];
        end
end

end

