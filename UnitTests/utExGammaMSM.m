classdef utExGammaMSM < utContinuous;
    
    properties (ClassSetupParameter)
        % Avoid similar rates for G and E, as this produces nearly singular info matrices in estimation.
        parmMuG    = struct( 'p200',200, 'p501',501 , 'p333',333, 'p170',170);
        parmSigmaG = struct( 'p80' , 80, 'p50' , 50 , 'p92' , 92, 'p30' , 30);
        parmMuE    = struct( 'p16' , 16, 'p72' , 72 , 'p111',111, 'p68' , 68);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExGammaMSM(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmMuG,parmSigmaG,parmMuE)
            % Computations specific to the ExGammaMSM distribution.
            testCase.Dist = ExGammaMSM(parmMuG,parmSigmaG,parmMuE);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExGammaMSM


