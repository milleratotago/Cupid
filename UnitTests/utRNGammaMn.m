classdef utRNGammaMn < utContinuous;
    
    properties (ClassSetupParameter)
        parmShape = struct( 'p3'  ,  3 , 'p10',10 , 'p1', 1, 'p21',21);
        parmMu    = struct( 'p150',150 , 'p9' , 9 , 'p2', 2, 'p53',53);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utRNGammaMn(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmShape,parmMu)
            % Computations specific to the RNGammaMn distribution.
            testCase.Dist = RNGammaMn(parmShape,parmMu);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utRNGammaMn


