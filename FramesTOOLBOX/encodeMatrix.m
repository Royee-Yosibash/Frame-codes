function [Mat_coded_partioned] = encodeMatrix(Mat, EncodingMat, rowOrCol)
%ENCODEMATRIX Summary of this function goes here
%
%
% Input
% -----------------------------
% (1) mat                   -   A matrix of arbitrary size to encode.
% (2) EncodingMat           -   The code generator matrix that should
%                               encode the 'mat' variable. The number of
%                               columns (m) defines the number of 
%                               partitions in 'mat' and the number of the  
%                               rows (n) defines the number of martirces at 
%                               the output.
% (3) rowOrCol              -   A string of either 'col' or 'row' that
%                               should decide if the partition should be
%                               done by columns of rows.
%
% Output:
% -----------------------------
% (1) Mat_coded_partioned	-	A cell vector of length n. In each cell one
%                               of the encoded partitions of the matrix in 
%                               'mat'. 
%
% Created by Royee Yosibash.
% Royeeyosibash@hotmail.com
% September 2021

[n,m] = size(EncodingMat);
Mat_partioned = cell(1,m);
Mat_coded_partioned = cell(1,n);
siz = size(Mat);

%% Partion the input matrices
switch rowOrCol
    case 'col'
        sizePartioned = siz(2)/m;
        if round(sizePartioned) ~= sizePartioned
            error('Matrix does not fit equal partioning');
        end
        
        for i = 1:m
            startPartionIndex = 1 + (i-1)*sizePartioned;  EndPartionIndex = i*sizePartioned;
            Mat_partioned{i} = Mat(:, (startPartionIndex:EndPartionIndex));
        end
    case 'row'
        sizePartioned = siz(1)/m;
        if round(sizePartioned) ~= sizePartioned
            error('Matrix does not fit equal partioning');
        end
        
        for i = 1:m
            startPartionIndex = 1 + (i-1)*sizePartioned;  EndPartionIndex = i*sizePartioned;
            Mat_partioned{i} = Mat((startPartionIndex:EndPartionIndex), :);
        end
end


%% Encode the Matrix
for i = 1:n
    EncodingRow = EncodingMat(i,:);
    Mat_coded_partioned{i}= 0;
    for j = 1:m
        Mat_coded_partioned{i} = Mat_coded_partioned{i} + EncodingRow(j)*Mat_partioned{j};
    end
end

end

