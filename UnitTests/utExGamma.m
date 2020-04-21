classdef utExGamma < utContinuous;
    
    properties (ClassSetupParameter)
        % Avoid similar rates for G and E, as this produces nearly singular info matrices in estimation.
        parmN     = struct( 'p1_5',1.5   , 'p1',1   , 'p3',3, 'p10',10);
        parmRateG = struct( 'p_008',.008 , 'p_1',.1 , 'p2',2, 'p10',10);
        parmRateE = struct( 'p_016',.016 , 'p_2',.2 , 'p1',1, 'p6',6);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExGamma(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmRateG,parmRateE)
            % Computations specific to the ExGamma distribution.
            testCase.Dist = ExGamma(parmN,parmRateG,parmRateE);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExGamma


