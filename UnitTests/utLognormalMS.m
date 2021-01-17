classdef utLognormalMS < utContinuous;
    
    properties (ClassSetupParameter)
        parmpostmu    = struct( 'p2',2 , 'p10',10 , 'p50',50 , 'p200',200 );
        parmpostsigma = struct( 'p1',1 , 'p2',2   , 'p5',5   , 'p50',50   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utLognormalMS(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmpostmu,parmpostsigma)
            % Computations specific to the distribution.
            testCase.Dist = LognormalMS(parmpostmu,parmpostsigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,2000);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmpostmu;
            testCase.Expected.SD = parmpostsigma;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.MLParmTolSE = 1.25;   % ML parameter estimation is not great
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.05;
            % Numerical problems recovering Kurtosis from CenMoments, especially with larger parmpostsigma.
            testCase.KurtRelTol = max( testCase.KurtRelTol, 0.03*parmpostsigma^2);

            testCase.SkipMGFs = true;
            testCase.SkipEstMoment = true;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utLognormalMS


