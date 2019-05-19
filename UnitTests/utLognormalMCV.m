classdef utLognormalMCV < utContinuous;
    
    properties (ClassSetupParameter)
        parmpostmu = struct( 'p2',2  , 'p10',10, 'p50',50, 'p200',200 );
        parmpostcv = struct( 'p_1',.1, 'p_2',.2, 'p_5',.5, 'p_3',.3   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utLognormalMCV(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmpostmu,parmpostcv)
            % Computations specific to the distribution.
            testCase.Dist = LognormalMCV(parmpostmu,parmpostcv);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,2000);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmpostmu;
            testCase.Expected.SD = parmpostcv*parmpostmu;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.MLParmTolSE = 1.25;   % ML parameter estimation is not great
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.05;
            % Numerical problems recovering Kurtosis from CenMoments, especially with larger parmpostcv.
            testCase.KurtRelTol = max( testCase.KurtRelTol, 0.03*parmpostcv^2);

            testCase.SkipMGFs = true;
            testCase.SkipMomentEst = true;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utLognormalMCV


