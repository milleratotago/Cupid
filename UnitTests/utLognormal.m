classdef utLognormal < utContinuous;
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'n2',-2    , 'n1',-1  , 'p0',0 , 'p1',1     ); % , 'p5',5     );  MOMENTS INACCURATE
        parmsigma = struct( 'p_01',.01 , 'p_1',.1 , 'p1',1 , 'p1_2',1.2 ); % , 'p2_1',2.1 );  FOR LARGER VALUES
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utLognormal(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmsigma)
            % Computations specific to the Lognormal distribution.
            testCase.Dist = Lognormal(parmmu,parmsigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.MLParmTolSE = 1.25;   % ML parameter estimation is not great
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.05;
            % Numerical problems recovering Kurtosis from CenMoments, especially with larger parmsigma.
            testCase.KurtRelTol = max( testCase.KurtRelTol, 0.03*parmsigma^2);

            testCase.SkipMGFs = true;
            testCase.SkipMomentEst = true;    % NEWJEFF: Not sure why this is here.  Under some situations it is possible to estimate with the method of moments.

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utLognormal


