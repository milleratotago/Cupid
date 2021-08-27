classdef utLogitTransDisc < utDiscrete
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utLogitTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LogitTrans distribution.
            % Basis distributions must be defined over (0-1) range.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LogitTrans(LinearTrans(Binomial(53,.44),1/154,.25));
                    testCase.EstParmCodes = 'frff';
                case 2
                    testCase.Dist = LogitTrans(LinearTrans(BinomialMixed([(1:10)/15]),0.07,0.05));
                    testCase.EstParmCodes = 'ffffffffffr';  % P's cannot be varied but they are counted as parameters
                    testCase.SkipEstAll = true;
                case 3
                    testCase.Dist = LogitTrans(MultTrans(Geometric(0.1),1e-8));
                    testCase.EstParmCodes = 'rf';
                case 4
                    p = Poisson(23.5);
                    testCase.Dist = LogitTrans(LinearTrans(p,1/(1.5*p.UpperBound),.01));
                    testCase.EstParmCodes = 'frf';
                case 5
                    testCase.Dist = LogitTrans(LinearTrans(UniformInt(0,100),.005,.01));
                    testCase.EstParmCodes = 'ffrf';
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
    
end  % utLogitTransDisc


