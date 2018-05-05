classdef utUniformInt < utDiscrete;
    
    properties (ClassSetupParameter)
        parmLo = struct('p2',2,   'p4',4   , 'p20',20, 'p100',100 , 'p500',500 );
        parmHi = struct('p10',10, 'p83',83 , 'p80',80, 'p550',550 , 'p730',730 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utUniformInt(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmLo,parmHi)
            % Computations specific to the normal distribution.
            testCase.Dist = UniformInt(parmLo,parmHi);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.LowerBound:testCase.Dist.UpperBound;
            
            % Set up some X values for which MLE should return the true parameters:
            npoints = 500;
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = (parmLo+parmHi)/2;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            % testCase.KurtRelTol = 0.01;  % Kurtosis not so accurate
            
            testCase.SkipAllEst = true;  % error when fminsearch suggests low/high crossover

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utUniformInt


