classdef utExtrVal1 < utContinuous;
    
    properties (ClassSetupParameter)
        parmalpha = struct( 'p_5',.5 , 'p1',1   , 'p3',3, 'p10',10 , 'p500',500 );
        parmbeta  = struct( 'p5',5   , 'p_1',.1 , 'p2',2, 'p10',10 , 'p100',100 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExtrVal1(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmalpha,parmbeta)
            % Computations specific to the ExtrVal1 distribution.
            testCase.Dist = ExtrVal1(parmalpha,parmbeta);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExtrVal1


