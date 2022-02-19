classdef FrameParameters < handle
    %FRAMEPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    %
    % Created by Royee Yosibash.
    % Royeeyosibash@hotmail.com
    % September 2021

    properties
        Type                % A string. This property indicates the type of frame that is described in this class
        SubType             % If the frame type has a subtype, this is the property were this is indicated
        M                   % A vector of positive integer values. Determines the number of rows in the frame. 
        N                   % A single value, a vector the same size as M. Determines the number of columns in the frame. 
        Gamma               % The ratios between the sets M to N
        normDim             % A string. The dimension that the frame should be normilized on. See the constructor function's documentation for supported inputs.
        Matrices            % This property stores previously created matrices. Used mainly to spare the runtime of recreating the same large matrices.
        isDeterministic     % A bolean that indicates if the the frames indicated above have only one possible value
        numPossibilities	% Counts the number different frames that fit the descriptions of every set of values. If the value is too large, the value is 'Inf'
    end
    
    methods
        %% Constructor fucntion
        function obj = FrameParameters(Type, SubType, M, N, Gamma, normDim)
            %FRAMEPARAMETERS Construct an instance of this class. 
            % Note: It is recommended to have the GUI construct members of 
            % this class. Use with caution if done so manualy.
            %
            %	Input
            % ------------------------
            % All inputs for this functions are the properties for this 
            % class type.
            % (1) Type          -   A string. See the property of the same
            %                       name for further explanation.
            %
            % (2) SubType       -   A string. See the property of the same
            %                       name for further explanation.
            %
            % (3) M             -   A vector of positive integer values. 
            %                       See the property of the same name for 
            %                       further explanation.
            %
            % (3) N             -   A vector of positive integer values. 
            %                       See the property of the same name for 
            %                       further explanation.
            % 
            % (4) Gamma         -   A vector of positive integer values. 
            %                       See the property of the same name for 
            %                       further explanation.
            %
            % (5) normDim       -   One of the following strings: 
            %                       a. 'None'
            %                       b. 'Row'
            %                       c. 'Column'
            %
            %   Output
            % ------------------------
            % (1) obj           -	A member of this class with the inputed 
            %                       parameters.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            availableNormDimStr = {'None'; 'Row'; 'Column'};
            FrameConstructorParseInput(Type, SubType, M, N, Gamma, normDim, availableNormDimStr);
            
            obj.Type = Type;
            obj.SubType = SubType;
            obj.M = M;
            obj.N = N;
            obj.Gamma = Gamma;
            obj.resetFrames();
            
        end

        
        %% 'Public' functions
        function [results] = gatherStatisticsfromAllFrames(obj, P, numTests, performanceMeasure , alpha, varargin)
            % This function returns a set of vectors each corresponds to
            % each of the frames in the obj. Each vector shows all
            % calculated performance measures as inputed by the user.
            %	Input
            % ------------------------
            % (1) obj                   -   A member of this class.
            %
            % (2) P                     -   A double variable in the range
            %                               of [0,1]. Describes the
            %                               probability that a vector in 
            %                               the over complete basis that 
            %                               constructs the frame is 
            %                               retained and not erased.
            %
            % (3) numTests              -   A single integer. The number of
            %                               times to run the tests to
            %                               gather statistics.
            %
            % (4) performanceMeasure	-	The required performance to
            %                               evaluate the frame the
            %                               available measures are:
            %                               a.'Condition Number'
            %                               b. 'Alpha Truncate'
            %                               c. 'Log Average'
            %                               d. 'Eigenvalue Offset +1'
            % 
            % (5) alpha                 -   A double variable in the range
            %                               of [0,1]. Describes the alpha
            %                               precentage to remove if 'Alpha
            %                               Truncate' is selected as the
            %                               performance measure. Note that
            %                               alpha is the precentage that is
            %                               removed, not retained.
            %
            % (6) varargin              -   Currently no additional inputs
            %                               are supported.
            %
            %   Output
            % ------------------------
            % (1) results               -	A vector of numeric values.    
            %                               Each element corresponding to a   
            %                               described frame and its value  
            %                               is the result to the
            %                               performance measure required.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.resetFrames();
            num = numel(obj.N); results = cell(1,num);
            beta = obj.Gamma/P; % m/k
            for idx = 1:num
                results{idx} = obj.gatherStatisticsfromSingleFrame(idx, P, numTests, performanceMeasure, alpha);
            end
        end

        
        function [gammaInv, result] = getNoiseAmplificationPerGamma(obj, P, numTests, PerformanceMeasure , alpha, varargin)
            % This function will be depriciated in future versions. Please
            % use gatherStatisticsfromAllFrames to recieve the results to
            % any performance query. the gammaInv vector can be calculated
            % by using gammaInv = obj.Gamma.^-1;
            %
            %	Input
            % ------------------------
            % (1) obj                   -   A member of this class.
            %
            % (2) P                     -   A double variable in the range
            %                               of [0,1]. Describes the
            %                               probability that a vector in 
            %                               the over complete basis that 
            %                               constructs the frame is 
            %                               retained and not erased.
            %
            % (3) numTests              -   A single integer. The number of
            %                               times to run the tests to
            %                               gather statistics.
            %
            % (4) performanceMeasure	-	The required performance to
            %                               evaluate the frame the
            %                               available measures are:
            %                               a. 'Condition Number'
            %                               b. 'Alpha Truncate'
            %                               c. 'Log Average'
            %                               d. 'Eigenvalue Offset +1'
            % 
            % (5) alpha                 -   A double variable in the range
            %                               of [0,1]. Describes the alpha
            %                               precentage to remove if 'Alpha
            %                               Truncate' is selected as the
            %                               performance measure. Note that
            %                               alpha is the precentage that is
            %                               removed, not retained.
            %
            % (6) varargin              -   Currently no additional inputs
            %                               are supported.
            %
            %   Output
            % ------------------------
            %
            % (1) gammaInv              -   A vector of numeric values that
            %                               represent the inversed values 
            %                               of gamma.
            %
            % (2) results               -	A vector of numeric values.    
            %                               Each element corresponding to a   
            %                               described frame and its value  
            %                               is the result to the
            %                               performance measure required.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.resetFrames();
            gammaInv = obj.Gamma.^-1;
            result = zeros(size(gammaInv));
            for idx = 1:numel(result)
                results = obj.gatherStatisticsfromSingleFrame(idx, P, numTests, PerformanceMeasure, alpha);
                result(idx) = min(results);
            end
        end
        
        
        %% 'Private' functions
        function [results] = gatherStatisticsFrameByFrame(obj, P, numTests, PreformanceMeasure , alpha, varargin)
            % This function returns a set of vectors each corresponds to
            % each of the frames in the obj. Each vector shows all
            % available noise amplification for that matrix permutation.
            obj.resetFrames();
            num = numel(obj.N); results = cell(1,num);
            for idx = 1:num
                results{idx} = obj.gatherStatisticsfromSingleFrame(idx, P, numTests, PreformanceMeasure, alpha);
            end
        end      
                
        
        function resetFrames(obj)
            % Cleanly resets the properties in this class that are used 
            % before a new run. 
            obj.Matrices = cell(size(obj.N));
            obj.numPossibilities = zeros(size(obj.N));
        end
        
        
        function generateSingleMatrix(obj, idx, varargin)
            % Creates a single matrix that fits the frame parameter in the
            % 'idx' index of the parameters in this class.
            [Code, obj.isDeterministic(idx), obj.numPossibilities(idx)] = ...
                getCode(obj.N(idx), obj.M(idx), obj.Type, obj.SubType, varargin{:});
            
            switch obj.normDim
                case 'Column'
                    for i = 1:size(Code,2)
                        Code(:,i) = Code(:,i)./norm(Code(:,i));
                    end
                    
                case 'Row'
                    for i = 1:size(Code,1)
                        Code(i,:) = Code(i,:)./norm(Code(i,:));
                    end
                case 'None'
            end
            CodeTr = transpose(Code);
            
            decoder = inv(ctranspose(CodeTr)*CodeTr) * ctranspose(CodeTr);
            
%             switch obj.normDim
%                 case 'Row'
%                     for i = 1:size(decoder,1)
%                         decoder(i,:) = decoder(i,:)./norm(decoder(i,:));
%                     end
%                 case 'Column'
%                     for i = 1:size(decoder,2)
%                         decoder(:,i) = decoder(:,i)./norm(decoder(:,i));
%                     end
%                 case 'None'
%             end
            
            obj.Matrices(idx) = {decoder};

        end
        
        function generateAllMatrixs(obj)
            % Generates all matrices for the the frame type described in
            % by this class.
            obj.isDeterministic = zeros(size(obj.N));
            obj.numPossibilities = obj.isDeterministic;
            obj.Matrices = cell(size(obj.N));
            for i = 1:numel(obj.Matrices)
                obj.generateSingleMatrix(i);
            end
        end
        
        
        function [H, isDeterministic] = getMatrix(obj, idx)
            % A query that gets the pre-existing matrix if it exists and
            % valid. If none exists or is valid, a valid matrix is created
            % and stored for further uses.

            % Check if the matrices have been initialized before
            if isempty(obj.Matrices)
                obj.generateAllMatrixs();
                % TODO: return the right matrix now because it was freshly
                %       created
            elseif isempty(cell2mat(obj.Matrices(1)))
                obj.generateAllMatrixs();
                % TODO: Same here
            end
            
            isDeterministic = obj.isDeterministic(idx);
            if ~isDeterministic
                [H, isDeterministic] = getCode(obj.N(idx), obj.M(idx), obj.normDim, obj.Type, obj.SubType); % TODO: Is this right? I think I removed the normDim
            else
                H = obj.Matrices{idx};
            end
        end
        
        
        function results = getAllPerformanceMeasureValues(obj, idx, beta, numTests, PerformanceMeasure, alpha)
            %	Input
            % ------------------------
            % (1) obj                   -   A member of this class.
            %
            % (2) idx                   -   An integer greater than 0. The 
            %                               index of the parameters of the 
            %                               frame that is under test.
            %
            % (3) beta                  -   A double variable that is
            %                               greater than 0. Describes the
            %                               ratio between the length of the
            %                               vectors in the over complete 
            %                               basis that constructs the frame 
            %                               to the number of vectors in
            %                               that basis after the erasures
            %                               of some of the vectors.
            %
            % (4) numTests              -   A single integer. The number of
            %                               times to run the tests to
            %                               gather statistics.
            %
            % (5) performanceMeasure	-	The required performance to
            %                               evaluate the frame the
            %                               available measures are:
            %                               a. 'Condition Number'
            %                               b. 'Alpha Truncate'
            %                               c. 'Log Average'
            %                               d. 'Eigenvalue Offset +1'
            % 
            % (6) alpha                 -   A double variable in the range
            %                               of [0,1]. Describes the alpha
            %                               precentage to remove if 'Alpha
            %                               Truncate' is selected as the
            %                               performance measure. Note that
            %                               alpha is the precentage that is
            %                               removed, not retained.
            %
            %	Output
            % --------------------------------------------------------
            %
            % (1) results	-	The performance measure values averaged 
            %                   over 'numTests' tests for the specified 
            %                   frame parameter at index 'idx'.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            ratio = zeros(size(obj.N));
            switch PerformanceMeasure
                case 'None'
                    [~,NoiseAmpValues] = collectMatrixEigenStatistics(obj, idx, beta, numTests, 0); % In natural units
                    results = NoiseAmpValues;
                case 'Condition Number'
                    [~,~, conditionNumbers] = collectMatrixEigenStatistics(obj, idx, beta, numTests, 0); % In natural units
                    results = conditionNumbers;
                case 'Alpha Truncate'
                    [~,NoiseAmpValues] = collectMatrixEigenStatistics(obj, idx, beta, numTests, 0); % In natural units
                    NoiseAmpValues = sort(NoiseAmpValues);
                    num2Del = floor(numel(NoiseAmpValues) * (alpha/2));
                    idx2Keep = (num2Del+1):(numel(NoiseAmpValues)-num2Del);
                    idx2Del = [(1:num2Del), ((numel(NoiseAmpValues)-num2Del+1):numel(NoiseAmpValues))];
                    figure(); hold on;
                    plot(idx2Del,10*log10(NoiseAmpValues(idx2Del)), 'ro');
                    plot(idx2Keep,10*log10(NoiseAmpValues(idx2Keep)), 'bo');
                    ratio(idx) = mean(NoiseAmpValues)/mean(NoiseAmpValues(idx2Keep));
                    results = NoiseAmpValues(idx2Keep);
                case 'Log Average'
                    [~,NoiseAmpValues] = collectMatrixEigenStatistics(obj, idx, beta, numTests, 0); % In natural units
                    results = 10*log10(NoiseAmpValues);
                case 'Eigenvalue Offset +1'
                    [~,NoiseAmpValues] = collectMatrixEigenStatistics(obj, idx, beta, numTests, 1); % In natural units
                    results = NoiseAmpValues;
                otherwise
                    error('getAllPerformanceMeasureValues: The Performance Measure specifeid is not implemented');
            end
        end
        
        
        function [meanResults] = gatherStatisticsfromSingleFrame(obj, idx, P, numTests, PerformanceMeasure, alpha)
            % Sets up and runs the required test for the user's performance 
            % measure multiple times for each possible frame described in
            % the class in the index 'idx'.
            maxPosibilites = 200;
            nPosibililties = min(maxPosibilites, obj.numPossibilities(idx));
            beta = obj.Gamma(idx)/P;
            obj.generateSingleMatrix(idx);  
            prevState = obj.isDeterministic(idx);
            
            meanResults = zeros(1,nPosibililties);
            for jPossibilty = 1:nPosibililties
                if mod(jPossibilty, 10) == 0
                    dispStr = ['Frame #',  num2str(idx),': Running variant ' , num2str(jPossibilty), ' of ', num2str(nPosibililties)];
                    disp(dispStr);
                end
                obj.generateSingleMatrix(idx, 'permIdx', jPossibilty);
                obj.isDeterministic(idx) = true;
                
                
                AllResults = getAllPerformanceMeasureValues(obj, idx, beta, numTests, PerformanceMeasure, alpha);
                meanResults(jPossibilty) = mean(AllResults);
            end
            
            
            if strcmpi(PerformanceMeasure,'Condition Number')
                % Doesnt need to be in dB
                bestResultStr = num2str(min(meanResults));
                disp(['Min condition number for idx ', num2str(idx) ,'# is ' , bestResultStr]);
            else
                if strcmpi(PerformanceMeasure,'Log Average')
                    % Already dB or doesnt need to be in dB
                else
                    meanResults  = 10*log10(meanResults); % to dB
                    if meanResults(jPossibilty) < 0 
                        error('stop here');
                    end
                end
                bestResultStr = num2str(min(meanResults));
                disp(['Min value for idx ', num2str(idx) ,'# is ' , bestResultStr, 'dB']);
            end
            
            obj.isDeterministic(idx) = prevState;
        end
        
    end
end


%% Private functions
function FrameConstructorParseInput(Type, SubType, M, N, Gamma, normDim, availableNormDimStr)
    p = inputParser;
    validationTypeFcn = @(x) isstring(x);
    validationSubTypeFcn = @(x) isstring(x);
    validationMFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validationNFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0); 
    validationGammaFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 1);
    validationnormDimFcn = @(x) isstring(x) && any(strcmpi(x, availableNormDimStr));
    
    addRequired(p,'Type',validationTypeFcn);
    addRequired(p,'SubType',validationSubTypeFcn);
    addRequired(p,'M',validationMFcn);
    addRequired(p,'N',validationNFcn);
    addRequired(p,'Gamma',validationGammaFcn);
    addRequired(p,'normDim',validationnormDimFcn);
    
    parse(p, Type, SubType, M, N, Gamma, normDim)
end


function [eigenvaluesVect, NoiseAmpValues, conditionNumbers] = collectMatrixEigenStatistics(Frame, idx, beta, numTests, eigvalueOffset)
%	Input
% ------------------------
% (1) obj                   -   A member of this class.
%
% (2) idx                   -   An integer greater than 0. The index of the 
%                               parameters of the frame that is under test.
%
% (3) beta                  -   A double variable that is
%                               greater than 0. Describes the
%                               ratio between the length of the
%                               vectors in the over complete 
%                               basis that constructs the frame 
%                               to the number of vectors in
%                               that basis after the erasures
%                               of some of the vectors.
%
% (4) numTests              -   A single integer. The number of
%                               times to run the tests to
%                               gather statistics.
%
% (5) eigvalueOffset        -   A single double variable. If the eigenvalue
%                               need to have an offset, it is added using 
%                               this variable. Note that the offset is
%                               calculated before the amplification or
%                               condition number, so these outputs might
%                               not be calculated in the way the were
%                               supposed to.
%
%   Output
% --------------------------------------------------------
%
%	(1) eigenvalVect        -	All Gram matrix eigenvalues from the  
%                               multiple random matrices.
%
%	(2) NoiseAmpValues      -	All noise amplification values calulated  
%                               for each test.
%
%   (3) conditionNumbers	-   All condition numbers calulated from the 
%                               frame entered.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = Frame.M(idx); n = Frame.N(idx); k = round(m/beta);
numEigen = m;    
eigenMat = inf * ones(numTests, numEigen);
NoiseAmpValues = zeros(numTests, 1);    conditionNumbers = NoiseAmpValues;
waitbarMsg = ['Running ' , num2str(numTests) ,' tests...'];  hWaitbar = waitbar(0, waitbarMsg);
for jTest = 1:numTests
    waitbar((jTest/numTests), hWaitbar, waitbarMsg);
    H = Frame.getMatrix(idx);
    idx = sort(randperm(n,k)); % choose random k of n
    if size(H,1) > size(H,2)
        H_idx = H(idx,:);
    else
        H_idx = H(:, idx);
    end
    
    [eigenValues, conditionNumbers(jTest)] = getGramMatrixEigenvalues(H_idx); % these are the eigenvalues of the gram matrix
    if any(eigenValues < 0 )
        NoiseAmpValues(jTest) = intmax; 
    else
        eigenMat(jTest, :) = eigenValues + eigvalueOffset;
        NoiseAmpValues(jTest) = sum(1./eigenMat(jTest, :))*sum(eigenMat(jTest, :))/(m^2);
    end
    if NoiseAmpValues(jTest)< 0 
        error('stop here')
    end
    
%     conditionNumbers(jTest) = abs(max(eigenMat(jTest,:))/min(eigenMat(jTest, :)));  % this should be taken as the root, as condition number of the gram matrix is ^2
end
close(hWaitbar);

conditionNumbers(conditionNumbers < 0) = inf;


eigenvaluesVect = reshape(eigenMat, [1, numel(eigenMat)]);
end

