classdef utList < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end

    methods
        
        function testCase=utList(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to this distribution.
            testCase.ThisCase  = parmCase;
            testCase.EstParmCodes = '';
            testCase.SkipAllEst = true;
            switch parmCase
                case 1
                    Xs = [2 4 8 100];
                    Ps = [2 4 8 100];
                case 2
                    Xs = 0.1:0.3:2.5;
                    Ps = Xs.^2;
                case 3
                    Xs = [-10 8 30 300];
                    Ps = [.2 .1 .4 .3];
            end
            Ps = Ps / sum(Ps);
            testCase.Dist = List(Xs,Ps);
            testCase.xvalues = testCase.Dist.DiscreteX;
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utList


