classdef utProductDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utProductDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Product distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = Product(UniformInt(1,5),UniformInt(10,20));
                    % testCase.EstParmCodes = 'rfff';
                    testCase.SkipEstAll = true;
                case 2
                    testCase.Dist = Product(Poisson(2.3),Binomial(10,.5));
                    testCase.EstParmCodes = 'ffr';
                case 3
                    testCase.Dist = Product(Geometric(.1),UniformInt(1,4));
                    testCase.EstParmCodes = 'rff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utProductDisc


