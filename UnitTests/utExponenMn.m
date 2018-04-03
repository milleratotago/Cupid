classdef utExponenMn < utContinuous;
    
    properties (ClassSetupParameter)
        parmmean  = struct( 'p205',205 , 'p97',97 , 'p11_3',11.3 , 'p1_7',1.7 , 'p_13',.13 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExponenMn(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmean)
            % Computations specific to the ExponenMn distribution.
            testCase.Dist = ExponenMn(parmmean);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,200);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmmean;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);  % Not very accurate
            testCase.MLParmTolSE = 0.25;   % ML parameter estimation is not great
            testCase.KurtRelTol = 0.004;
            testCase.CenMomentAbsTol(2) = max( testCase.CenMomentAbsTol(2), .00005*parmmean);  % Numerical problems at upper tail if rate small.
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExponenMn


