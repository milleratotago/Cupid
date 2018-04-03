classdef utRNGamma < utContinuous;
    
    properties (ClassSetupParameter)
        parmN    = struct( 'p_5',.5     , 'p1',1   , 'p3',3, 'p10',10);
        parmRate = struct( 'p_005',.005 , 'p_1',.1 , 'p2',2, 'p10',10);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utRNGamma(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmRate)
            % Computations specific to the RNGamma distribution.
            testCase.Dist = RNGamma(parmN,parmRate);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utRNGamma


