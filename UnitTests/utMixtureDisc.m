classdef utMixtureDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utMixtureDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Mixture distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = Mixture(0.3,UniformInt(1,5),UniformInt(10,20));
                    % testCase.EstParmCodes = 'rfff';
                    testCase.SkipEstAll = true;
                case 2
                    testCase.Dist = Mixture(0.7,Poisson(2.3),0.1,Binomial(10,.5),UniformInt(10,20));
                    testCase.EstParmCodes = 'frfffff';
                case 3
                    testCase.Dist = Mixture(0.9,Geometric(.1),UniformInt(1,4));
                    testCase.EstParmCodes = 'frff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utMixtureDisc


