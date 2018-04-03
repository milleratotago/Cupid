classdef utGamma < utContinuous;
    
    properties (ClassSetupParameter)
        parmN    = struct( 'p1',1     , 'p2',2   , 'p4',4, 'p10',10);
        parmRate = struct( 'p_005',.005 , 'p_1',.1 , 'p2',2, 'p10',10);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utGamma(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmRate)
            % Computations specific to the Gamma distribution.
            testCase.Dist = Gamma(parmN,parmRate);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Not numerically very accurate
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters.

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utGamma


