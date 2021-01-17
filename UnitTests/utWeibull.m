classdef utWeibull < utContinuous;
    
    properties (ClassSetupParameter)
        parmScale  = struct( 'p_9',.9   , 'p1',1 , 'p5',5     , 'p20',20   );
        parmPower  = struct( 'p_75',.75 , 'p1',1 , 'p1_5',1.5 , 'p1_2',1.2 );  % Numerical problems if power too far from 1
        parmOrigin = struct( 'p0',0     , 'p1',1 , 'p5',5     , 'n10',-10  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utWeibull(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmScale,parmPower,parmOrigin)
            % Computations specific to the Weibull distribution.
            testCase.Dist = Weibull(parmScale,parmPower,parmOrigin);
            testCase.EstParmCodes = 'rrf';
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,1000,10000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            % Parameter estimates are not very accurate.
            testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.10;
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.05;

            testCase.SkipEstMoment = true;  % Moment estimation is pretty hopeless for Weibull

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utWeibull


