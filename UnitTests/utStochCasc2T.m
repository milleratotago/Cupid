classdef utStochCasc2T < utContinuous
    
    properties (ClassSetupParameter)
        parmrate1  = struct( 'p_005',.005 , 'p_012',.012 , 'p_01',.01    , 'p1',1    );
        parmrate2  = struct( 'p_015',.015 , 'p_065',.065 , 'p_1',.1      , 'p2',2    );
        ParmK      = struct( 'p_2',2      , 'p_3' ,3     , 'p_10',10     , 'p_50',50 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utStochCasc2T(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmrate1,parmrate2,ParmK)
            % Computations specific to the StochCasc2T distribution.
            testCase.Dist = StochCasc2T(parmrate1,parmrate2,ParmK);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);  % Not very accurate1
            testCase.ParmEstAbsTol = 0.02;
            testCase.ParmEstRelTol = 0.02;
            testCase.MLParmTolSE = 1.5;
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utStochCasc2T


