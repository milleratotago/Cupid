classdef utLoglogistic < utContinuous;
    
    properties (ClassSetupParameter)
        parmscale = struct( 'p200',200 , 'p275',275  , 'p500',500 , 'p1000',1000 );
        parmshape = struct( 'p4',4     , 'p3',3      , 'p4_1',4.1 , 'p5_2',5.2   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utLoglogistic(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmscale,parmshape)
            % Computations specific to the Loglogistic distribution.
            testCase.Dist = Loglogistic(parmscale,parmshape);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

%           testCase.SkipMGFs = true;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utLoglogistic


