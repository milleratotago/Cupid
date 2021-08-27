classdef utVonMises < utContinuous
    
    properties (ClassSetupParameter)
        parmloc   = struct( 'p_1',.1   , 'p_85',.85  , 'p1_1',1.1 , 'p1_5',1.5   , 'p3_3',3.3   , 'p6_2',6.2 );
        parmConcentration = struct( 'p_01',.01 , 'p_1',.1    , 'p1',1     , 'p2',2       , 'p4',4     , 'p11',11   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utVonMises(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmloc,parmConcentration)
            % Computations specific to the VonMises distribution.
            testCase.Dist = VonMises(parmloc,parmConcentration);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utVonMises


