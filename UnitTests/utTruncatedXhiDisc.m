classdef utTruncatedXhiDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utTruncatedXhiDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the TruncatedXhi distribution.
            % Basis distribution values should not be too large
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = TruncatedXhi(Binomial(53,.84),55);
                    testCase.EstParmCodes = 'frf';
                case 2
                    testCase.Dist = TruncatedXhi(BinomialMixed([(1:10)/11]),5);
%                     testCase.EstParmCodes = 'ffffffffffr';  % P's cannot be varied but they are counted as parameters
                    testCase.SkipEstAll = true;
                case 3
                    testCase.Dist = TruncatedXhi(Geometric(0.7),5);
                    testCase.EstParmCodes = 'rf';
                case 4
                    testCase.Dist = TruncatedXhi(Poisson(23.5),35);
                    testCase.EstParmCodes = 'rf';
                case 5
                    testCase.Dist = TruncatedXhi(UniformInt(1,25),20);  % This case checked in ComparePDF
%                     testCase.EstParmCodes = 'fff';
                    testCase.SkipEstAll = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.XGrain = 10;  % Be lenient on accepting X values
            testCase.Dist.CDFGrain = 10;  % Be lenient on accepting X values
            
            testCase.SetupXs;  % xvalues at which PDF, CDF, etc should be evaluated & xMLE used to estimate MLE

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 5
                tempU = UniformInt(1,20);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end % methods (Test)
    
end  % utTruncatedXhiDisc


