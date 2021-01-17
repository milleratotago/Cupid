classdef utArcsinTransDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utArcsinTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the ArcsinTrans distribution.
            % Basis distributions must be defined over 0-1 range.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ArcsinTrans(MultTrans(Binomial(53,.44),1/54));
                    testCase.EstParmCodes = 'frf';
                case 2
                    testCase.Dist = ArcsinTrans(LinearTrans(BinomialMixed([(1:10)/11]),0.8,.1));
                    testCase.EstParmCodes = 'ffffffffffrf';  % P's cannot be varied but they are counted as parameters
                    testCase.SkipEstAll = true;
                case 3
                    testCase.Dist = ArcsinTrans(MultTrans(Geometric(0.1),.001));
                    testCase.EstParmCodes = 'fr';
                case 4
                    testCase.Dist = ArcsinTrans(MultTrans(Poisson(23.5),.005));
                    testCase.EstParmCodes = 'fr';
                case 5
                    testCase.Dist = ArcsinTrans(LinearTrans(UniformInt(0,100),1/200,.5));
                    testCase.EstParmCodes = 'fffr';
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
    
end  % utArcsinTransDisc


