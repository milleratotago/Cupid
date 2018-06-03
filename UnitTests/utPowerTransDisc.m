classdef utPowerTransDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utPowerTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the PowerTrans distribution.
            % Basis distribution values should not be too large
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = PowerTrans(Binomial(53,.84),.37);
                    testCase.EstParmCodes = 'frf';
                case 2
                    testCase.Dist = PowerTrans(BinomialMixed([(1:10)/11]),1);
%                     testCase.EstParmCodes = 'ffffffffffr';  % P's cannot be varied but they are counted as parameters
                    testCase.SkipAllEst = true;
                case 3
                    testCase.Dist = PowerTrans(Geometric(0.7),1.3);
                    testCase.EstParmCodes = 'rf';
                case 4
                    testCase.Dist = PowerTrans(Poisson(23.5),1.25);
                    testCase.EstParmCodes = 'rf';
                case 5
                    testCase.Dist = PowerTrans(UniformInt(1,25),1.52);
%                     testCase.EstParmCodes = 'ffr';
                    testCase.SkipAllEst = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.XGrain = 50;  % Be lenient on accepting X values
            
            testCase.SetupXs;  % xvalues at which PDF, CDF, etc should be evaluated & xMLE used to estimate MLE

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            testCase.ParmEstAbsTol = 0.05;
            testCase.ParmEstRelTol = 0.05;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utPowerTransDisc


