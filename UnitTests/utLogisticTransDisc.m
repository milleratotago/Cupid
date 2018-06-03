classdef utLogisticTransDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utLogisticTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LogisticTrans distribution.
            % Basis distributions must be defined over 0-1 range.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LogisticTrans(Binomial(9,.44));
                    testCase.EstParmCodes = 'fr';
                case 2
                    testCase.Dist = LogisticTrans(BinomialMixed([(1:10)/11]));
                    testCase.EstParmCodes = 'ffffffffff';  % P's cannot be varied but they are counted as parameters
                    testCase.SkipAllEst = true;
                case 3
                    testCase.Dist = LogisticTrans(Geometric(0.6));
                    testCase.EstParmCodes = 'r';
                case 4
                    testCase.Dist = LogisticTrans(Poisson(3.5));
                    testCase.EstParmCodes = 'r';
                case 5
                    testCase.Dist = LogisticTrans(UniformInt(0,10));
                    testCase.EstParmCodes = 'ff';
                    testCase.SkipAllEst = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.XGrain = 50;  % Be lenient on accepting X values
            
            testCase.SetupXs;  % xvalues at which PDF, CDF, etc should be evaluated & xMLE used to estimate MLE

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utLogisticTransDisc


