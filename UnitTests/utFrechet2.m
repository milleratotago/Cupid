classdef utFrechet2 < utContinuous
    
    properties (ClassSetupParameter)
        parmshape  = struct( 'p12_5',12.5 , 'p20_2',20.2, 'p8'   ,8    );
        parmscale  = struct( 'p_105',.105 , 'p10',10  , 'p2_25',2.25 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utFrechet2(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmshape,parmscale)
            % Computations specific to the Frechet2 distribution.
            testCase.Dist = Frechet2(parmshape,parmscale);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utFrechet2


