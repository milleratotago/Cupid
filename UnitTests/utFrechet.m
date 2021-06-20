classdef utFrechet < utContinuous;
    
    properties (ClassSetupParameter)
        parmshape  = struct( 'p2_5',2.5   , 'p1_2',1.2 , 'p8'   ,8    );
        parmscale  = struct( 'p_105',.105 , 'p10',10   , 'p2_25',2.25 );
        parmminval = struct( 'p0',0       , 'p2',2     , 'p5'   ,5    );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utFrechet(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmshape,parmscale,parmminval)
            % Computations specific to the Frechet distribution.
            testCase.Dist = Frechet(parmshape,parmscale,parmminval);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
           
%            testCase.MLParmTolSE = 2.5;   % ML parameter estimation is very bad
%            testCase.MGFMom2RelTol = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utFrechet


