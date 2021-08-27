classdef utNegativeBinomial < utDiscrete
    
    properties (ClassSetupParameter)
        parmN = struct('p2',2  , 'p4',4     , 'p20',20 );
        parmp = struct('p_1',.1, 'p_83',.83 , 'p_8',.8 );
    end
    % Note relaxed moment tolerances when using Poisson approximations
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utNegativeBinomial(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmp)
            % Computations specific to the distribution.
            testCase.Dist = NegativeBinomial(parmN,parmp);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.LowerBound:testCase.Dist.UpperBound;
            
            % Set up some X values for which MLE should return the true parameters:
            npoints = 500;
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utNegativeBinomial


