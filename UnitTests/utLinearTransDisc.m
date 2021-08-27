classdef utLinearTransDisc < utDiscrete
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utLinearTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LinearTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LinearTrans(Binomial(53,.44),-2.3,55);
                    testCase.EstParmCodes = 'frff';
                case 2
                    testCase.Dist = LinearTrans(BinomialMixed([(1:10)/11]),0.75,-5.5);
                    testCase.SkipEstAll = true;
%                     testCase.EstParmCodes = 'ffffffffffrr';  % P's cannot be varied but they are counted as parameters
                case 3
                    testCase.Dist = LinearTrans(Geometric(0.1),-2.3,5.5);
                    testCase.EstParmCodes = 'rff';  % Don't adjust both constant and mean since they can trade off.
                case 4
                    testCase.Dist = LinearTrans(Poisson(23.5),0.9,100);
                    testCase.EstParmCodes = 'rff';  % Don't adjust both constant and mean since they can trade off.
                case 5
                    testCase.Dist = LinearTrans(UniformInt(0,100),1.5,10);
                    testCase.SkipEstAll = true;
%                     testCase.EstParmCodes = 'ffrf';  % Constant to nearest integer because Xs are integers.
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.XGrain = 10;  % Be lenient on accepting X values
            
            testCase.SetupXs;  % xvalues at which PDF, CDF, etc should be evaluated & xMLE used to estimate MLE

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.015);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utLinearTransDisc


