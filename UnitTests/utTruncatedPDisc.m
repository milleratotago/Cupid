classdef utTruncatedPDisc < utDiscrete
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utTruncatedPDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the TruncatedP distribution.
            % Basis distribution values should not be too large
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = TruncatedP(Binomial(53,.84),.4,.9);
%                    testCase.EstParmCodes = 'frff';
                case 2
                    testCase.Dist = TruncatedP(BinomialMixed([(1:10)/11]),.2,.8);
%                     testCase.EstParmCodes = 'ffffffffffr';  % P's cannot be varied but they are counted as parameters
                case 3
                    testCase.Dist = TruncatedP(Geometric(0.7),.25,.90);
%                    testCase.EstParmCodes = 'rff';
                case 4
                    testCase.Dist = TruncatedP(Poisson(23.5),.15,.84);
                    testCase.EstParmCodes = 'rff';
                case 5
                    testCase.Dist = TruncatedP(UniformInt(1,25),.24,.80);  % Uniform assumed by ComparePDFs
%                     testCase.EstParmCodes = 'ffff';
                    testCase.SkipEstAll = true;
            end
            testCase.SkipEstAll = true;  % NWJEFF: Estimation may not really make sense because the P's aren't actually fixed
                                         % at the specified values, due to the graininess of the CDF function. In any case,
                                         % it does not recover the original parameter values.
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
                tempU = UniformInt(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end % methods (Test)
    
end  % utTruncatedPDisc


