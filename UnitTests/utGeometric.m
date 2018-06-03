classdef utGeometric < utDiscrete;
    
    properties (ClassSetupParameter)
        parmp = struct('p_1',0.1, 'p_33',.33 , 'p_55',.55 , 'p_73',.73, 'p_48',.48 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utGeometric(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmp)
            % Computations specific to the normal distribution.
            testCase.Dist = Geometric(parmp);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.SetupXs;
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Expect some problems due to truncation of discrete distribution

            if parmp>.6  % not enough probability in different bins to do Probit estimation
                testCase.SkipAllProbitEst = true;
            end
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utGeometric


