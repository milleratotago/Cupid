classdef utTruncatedX < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utTruncatedX(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the TruncatedX distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = TruncatedX(Uniform(100,110),102,107);  % Uniform assumed by ComparePDFs
                    testCase.SkipEstAll = true;  % The uniform bounds don't matter as long as they are outside the cutoffs.
                    testCase.Expected.Mean = (testCase.Dist.LowerBound + testCase.Dist.UpperBound)/2;
                case 2
                    testCase.Dist = TruncatedX(Normal(0,1),-1,1);
                    testCase.EstParmCodes = 'rrff';
                case 3
                    testCase.Dist = TruncatedX(Exponential(0.1),1,20);
                    testCase.EstParmCodes = 'rff';
                    % testCase.SkipEstML = true;  % The exponential produces a better ML at a different rate.
                case 4
                    testCase.Dist = TruncatedX(AddPosTrans(RNGammaMS(400,100),100),350,650);
                    testCase.EstParmCodes = 'rrrff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,50000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
                tempU = Uniform(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utTruncatedX


