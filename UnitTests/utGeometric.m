classdef utGeometric < utDiscrete;
    
    properties (ClassSetupParameter)
        parmp = struct('p_1',0.1, 'p_33',.33 , 'p_55',.55 , 'p_73',.73, 'p_84',.84 );
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

            testCase.xvalues = testCase.Dist.LowerBound:testCase.Dist.UpperBound;
            
            % Set up some X values for which MLE should return the true parameters:
            npoints = 500;
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            if parmp>.6  % not enough probability in different bins to do Probit estimation
                testCase.SkipAllProbitEst = true;
            end
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utGeometric


