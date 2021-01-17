classdef utInverseTransDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utInverseTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the InverseTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = InverseTrans(AddTrans(Binomial(35,.44),.25));
                    testCase.EstParmCodes = 'frf';
                case 2
                    testCase.Dist = InverseTrans(LinearTrans(BinomialMixed([(1:10)/15]),0.07,0.05));
%                     testCase.EstParmCodes = 'ffffffffffr';  % P's cannot be varied but they are counted as parameters
                    testCase.SkipEstAll = true;
                case 3
                    testCase.Dist = InverseTrans(AddTrans(Geometric(0.1),5));
                    testCase.EstParmCodes = 'rf';
                case 4
                    testCase.Dist = InverseTrans(AddTrans(Poisson(23.5),0.5));
                    testCase.EstParmCodes = 'rf';
                case 5
                    testCase.Dist = InverseTrans(UniformInt(1,10));
%                    testCase.EstParmCodes = 'fi';
                    testCase.SkipEstAll = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.SetupXs;  % xvalues at which PDF, CDF, etc should be evaluated & xMLE used to estimate MLE

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utInverseTransDisc


