classdef utNakaRush < utContinuous;
    
    properties (ClassSetupParameter)
        parmScale  = struct( 'p_1',.1 , 'p1',1 , 'p2',2 , 'p4',4 , 'p10',10 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utNakaRush(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmScale)
            % Computations specific to the NakaRush distribution.
            testCase.Dist = NakaRush(parmScale);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            
            SetupXs(testCase,41,500);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.MGFMom2RelTol = 0.05;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utNakaRush


