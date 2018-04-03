classdef utCauchy < utContinuous;
    
    properties (ClassSetupParameter)
        parmLoc   = struct( 'n100',-100 , 'p0',0 ,  'p10',10 , 'p100',100 );
        parmScale = struct( 'p0_5',0.5  , 'p5',5 ,  'p10',10 ,  'p25',25 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utCauchy(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmLoc,parmScale)
            % Computations specific to the Cauchy distribution.
            testCase.Dist = Cauchy(parmLoc,parmScale);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            % testCase.Dist.SearchOptions.TolFun = 1e-14;
            % testCase.Dist.SearchOptions.TolX = 1e-14;
            
            testCase.HighestMoment = 2;  % Need to compute at least 2 moments for EstMom,
            % but I will not attempt to check the computational accuracy
            % of variance, skewness, & kurtosis, which technically do not exist anyway.
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);  % Lots of numerical problems with Cauchy.
            testCase.CenMomentAbsTol(2) = max(testCase.CenMomentAbsTol(2), 0.010*(abs(testCase.Dist.Mean)+testCase.Dist.SD) );
            testCase.RawMomentAbsTol(3) = realmax;   % Do not check Moment2 accuracy
            testCase.MGFMom2RelTol = realmax;
            
            % Random numbers are unbounded.
            testCase.MinRand = -inf;
            testCase.MaxRand = inf;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utCauchy


