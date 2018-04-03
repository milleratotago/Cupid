classdef utQuantal < utContinuous;
    
    properties (ClassSetupParameter)
        parmthreshold  = struct( 'p1',1 , 'p2',2 , 'p4',4 , 'p8',8 , 'p16',16 , 'p32',32 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utQuantal(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmthreshold)
            % Computations specific to the Quantal distribution.
            testCase.Dist = Quantal(parmthreshold);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmthreshold;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);  % Not very accurate
%            testCase.MLParmTolSE = 0.25;   % ML parameter estimation is not great
%            testCase.KurtRelTol = 0.004;
%            testCase.CenMomentAbsTol(2) = max( testCase.CenMomentAbsTol(2), .00005/parmthreshold);  % Numerical problems at upper tail if rate small.
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters 
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utQuantal


