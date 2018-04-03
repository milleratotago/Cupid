classdef utRosin < utContinuous;
    
    properties (ClassSetupParameter)
        parmP80   = struct( 'p10',10   , 'p20',20 , 'p50',50 , 'p250',250 , 'p1250',1250 );
        parmpower = struct( 'p1_5',1.5 , 'p2_5',2.5 , 'p3',3   , 'p7',7     , 'p11',11   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utRosin(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmP80,parmpower)
            % Computations specific to the Rosin distribution.
            testCase.Dist = Rosin(parmP80,parmpower);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.RawMomentRelTol(4) = 0.005;
%             testCase.RawMomentRelTol(4) = 0.1;
%             testCase.MLParmTolSE = 0.2;
% 
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utRosin


