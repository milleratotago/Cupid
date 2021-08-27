classdef utPoisson < utDiscrete
    
    properties (ClassSetupParameter)
        parmmu = struct('p_8',0.8 , 'p5',5   , 'p250',25   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utPoisson(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu)
            % Computations specific to the Poisson distribution.
            testCase.Dist = Poisson(parmmu);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

%             SetupXs(testCase,100,2000);
            testCase.xvalues = testCase.Dist.XsToPlot;
            testCase.xMLE = parmmu*ones(20,1);  % Cheat to make sure mean of observations matches parmmu
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utPoisson


