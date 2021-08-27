classdef utrNoncentral < utContinuous
    
    properties (ClassSetupParameter)
        parmSampleSize = struct( 'p6',6       , 'p12',12 , 'p24',24 , 'p150',150 );
        parmTrueRho    = struct( 'p_005',.005 , 'p_3',.3 , 'p_5',.5 , 'p_1',.1   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utrNoncentral(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmSampleSize,parmTrueRho)
            % Computations specific to the tNoncentral distribution.
            testCase.Dist = rNoncentral(parmSampleSize,parmTrueRho);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Not numerically very accurate
            testCase.EstParmCodes = 'fr';  % Don't adjust the integer paramter
%             testCase.ParmEstAbsTol(2) = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utrNoncentral


