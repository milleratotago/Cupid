classdef utRNGammaMS < utContinuous
    
    properties (ClassSetupParameter)
        ParmMu = struct( 'p300',300 , 'p10',10 , 'p3',3, 'p210',210);
        ParmSD = struct( 'p150',150 , 'p9',9   , 'p2',2, 'p53',53);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utRNGammaMS(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,ParmMu,ParmSD)
            % Computations specific to the RNGammaMS distribution.
            testCase.Dist = RNGammaMS(ParmMu,ParmSD);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utRNGammaMS


