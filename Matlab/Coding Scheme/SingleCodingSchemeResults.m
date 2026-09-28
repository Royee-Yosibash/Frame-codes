classdef SingleCodingSchemeResults < handle
    %SINGLECODINGSCHEMERESULTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N
        M
        SNR
        resultsTable
        resultsTableLaTeX
    end
    
    methods
        function obj = SingleCodingSchemeResults(N, M, SNR, measuresRealName, measuresLatex)
            %SINGLECODINGSCHEMERESULTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.N                   = N;
            obj.M                   = M;
            obj.SNR                 = SNR;
            
            measurements = measuresRealName;
            obj.resultsTable        = table(measurements); %change to rownames?
            
            measurements = measuresLatex;
            obj.resultsTableLaTeX   = table(measurements);
        end
        
        function addTableCol(obj, k, codeType, ColValues)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Colname = num2str(k);
            Colname = [codeType,'_k_',Colname];
            Colname = strrep(Colname,' ','_');
            obj.resultsTable.(Colname) = ColValues;
            obj.resultsTableLaTeX.(Colname) = createStringsForLaTex(ColValues);
        end
        
        
        function [k,Vals] = getValuesForCode(obj, Measure, CodeType)
            codeStr = strrep(CodeType,' ','_');
%             indCodeType = strcmp(obj.CodeTypes , CodeType);
%             if ~any(indCodeType)
%                 error('No such codeType');
%             end
%             numInCodes = find(indCodeType);
          
            Table = obj.resultsTable;
            varNames = Table.Properties.VariableNames;
            idxCodes = contains(varNames,codeStr);
%             idxCodes = ~cellfun(@isempty,numInCodes);
            
            measurements = Table.measurements;
            indRow = strcmp(measurements,Measure);
            if ~any(indRow)
                error('No such measure');
            end
            
            newTable = Table(indRow,idxCodes);
            Vals = newTable{:,:}; Vals = cell2mat(Vals);
            
            k = newTable.Properties.VariableNames;
            for i=1:numel(k)
                temp = k{i};
                temp = strrep(temp, codeStr, '');
                temp = strrep(temp,'_k_', '');
                k{i} = str2num(temp);
            end
            
            k = cell2mat(k);
            
        end
        
        
        function haxes = plotByMeasure(obj, Measure, CodeType, in_axes)
            if ~isempty(in_axes)
                haxes = in_axes;
            else
                h = figure();
                haxes = axes(h);
            end
               
            [k,Vals] = obj.getValuesForCode(Measure, CodeType);
            
            linewidth = 2;  fontsize = 16;
            plot(k,Vals, '*-', 'DisplayName', CodeType, 'LineWidth', linewidth);
            xlabel('k', 'fontsize', fontsize);
            ylabel(Measure, 'fontsize', fontsize);
        end
        
        
        function writeToTXTfile(obj)
            %WRITETOTXT Summary of this function goes here
            %   Detailed explanation goes here
            fName = ['n = ', num2str(obj.N), ' m = ', num2str(obj.M), ' SNR = ', num2str(obj.SNR)];
            fileID = fopen([fName, '.txt'], 'w');

            TableToLatex = obj.resultsTableLaTeX;
            VarNames = TableToLatex.Properties.VariableNames;
            newLatexRow=[];
            for i= 1:numel(VarNames)
                newLatexRow = [newLatexRow, VarNames{i}, '\t'];
            end
            fprintf(fileID, [newLatexRow, '\r\n']);

            for i= 1:size(TableToLatex,1)
                TableRow = TableToLatex(i,:);
                newStr = cell2mat(TableRow{1,1});
                newLatexRow = strrep(newStr,'\','\\');

                for j = 2:size(TableRow,2)
                    newStr = cell2mat(TableRow{1,j});
                    newStr = strrep(newStr,'\','\\');        
                    newLatexRow = [newLatexRow, ' & $', newStr, '$'];
                end
                newLatexRow = [newLatexRow, ' \\\\\r\n'];
                fprintf(fileID, newLatexRow);
            end
            fclose(fileID);

        end
        
    end
end


function [out_strings] = createStringsForLaTex(ColValues)
%CREATESTRINGSFORLATEX Summary of this function goes here
%   Detailed explanation goes here

% temp = {num2str(MSE_mean); num2str(ForbeniusNormError_mean); num2str(meanCondNumber); num2str(minCondNumber); num2str(maxCondNumber)};

out_strings = cell(size(ColValues));

for i=1:numel(out_strings)
    out_strings{i} = num2str(ColValues{i});
end

for iTemp = 1:numel(out_strings)
    tempNum = str2num(out_strings{iTemp});
    flag = true; power = 0;
    while flag
        if abs(tempNum) < 10 && abs(tempNum) >= 1
            flag = false;
        elseif abs(tempNum) < 1
            tempNum = tempNum*10;
            power = power-1;
        else
            tempNum = tempNum/10;
            power = power+1;
        end
    end
    
    if power == 0
        powerString = [];
    elseif power == 1
        tempNum = tempNum*10;
        powerString = [];
    elseif power == -1
        tempNum = tempNum/10;
        powerString = [];
    else
        powerString = ['\cdot 10^{', num2str(power), '}'];
    end
    
    numStr = num2str(tempNum);
    if numel(numStr) > 4
        if tempNum > 0
            numStr = numStr(1:4);
        else
            numStr = numStr(1:5);
        end
        
    end
    totString = [numStr, powerString];
    out_strings{iTemp} = totString;
end

end


