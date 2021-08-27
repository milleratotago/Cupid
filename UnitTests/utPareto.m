classdef utPareto < utContinuous
    
    properties (ClassSetupParameter)
        parm1K = struct( 'p_5',.5  , 'p2',2 , 'p5',5   , 'p25',25 , 'p250',250 );
        parm2A = struct( 'p3',3    , 'p4',4 , 'p6',6   , 'p10',10 , 'p40',40   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utPareto(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1K,parm2A)
            % Computations specific to the Pareto distribution.
            testCase.Dist = Pareto(parm1K,parm2A);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,41,1500);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
%            testCase.RawMomentAbsTol(4) = 0.1;
%            testCase.RawMomentRelTol(4) = 0.1;
            testCase.KurtRelTol = 0.005;
            testCase.ParmEstRelTol = 0.01;
            testCase.MLParmTolSE = 1;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utPareto


