classdef utWald < utContinuous;
    
    
    properties (ClassSetupParameter)
        waldmu      = struct( 'p_20',.20 ,  'p_10',.10 ,  'p_40',.40  ,  'p_50',.50 );
        waldsigma   = struct( 'p1',1     ,  'p1b',1    ,  'p1c',1     ,  'p2',2     );
        waldbarrier = struct( 'p50',50   ,  'p35',35   ,  'p20',20    ,  'p45',45   );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utWald(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,waldmu,waldsigma,waldbarrier)
            % Computations specific to the Wald distribution.
            testCase.Dist = Wald(waldmu,waldsigma,waldbarrier);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            testCase.Dist.SearchOptions.MaxFunEvals = 20000;
            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,5000);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            % Parameter estimation is not very good
            testCase.ParmEstAbsTol = 0.06;
            testCase.ParmEstRelTol = 0.09;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utWald


