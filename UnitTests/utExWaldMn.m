classdef utExWaldMn < utContinuous;
    
    
    properties (ClassSetupParameter)
        Case          = struct( 'case1',1  , 'case2',2   , 'case3',3  , 'case4',4  );
        p1waldmu      = struct( 'p_20',.20 , 'p_20b',.20 , 'p_40',.40 , 'p_50',.50 );
        p2waldsigma   = struct( 'p1',1     , 'p1b',1     , 'p1c',1    , 'p2',2     );
        p3waldbarrier = struct( 'p50',50   , 'p50b',50   , 'p20',20   , 'p45',45   );
        p4ExMean      = struct( 'p5',5     , 'p10',10    , 'p15',15   , 'p30',30 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utExWaldMn(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,Case,p1waldmu,p2waldsigma,p3waldbarrier,p4ExMean)
            % Computations specific to the ExWaldMn distribution.
            testCase.ThisCase = Case;
            testCase.Dist = ExWaldMn(p1waldmu,p2waldsigma,p3waldbarrier,p4ExMean);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.EstParmCodes = 'rffr';
            testCase.Dist.UseSplinePDFOn(200);
            testCase.SkipMomentEst = true;  % Moments do not provide enough info to constrain parameters
            
            SetupXs(testCase,40,1000);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            % Parameter estimation is not very good
            testCase.ParmEstAbsTol = 0.06;
            testCase.ParmEstRelTol = 0.09;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utExWaldMn


