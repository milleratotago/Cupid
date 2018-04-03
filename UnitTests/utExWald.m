classdef utExWald < utContinuous;
    
    
    properties (ClassSetupParameter)
        Case        = struct( 'case1',1  ,  'case2',2   ,  'case3',3   ,  'case4',4  );
        waldmu      = struct( 'p_20',.20 ,  'p_20b',.20 ,  'p_40',.40  ,  'p_50',.50 );
        waldsigma   = struct( 'p1',1     ,  'p1b',1     ,  'p1c',1     ,  'p2',2     );
        waldbarrier = struct( 'p50',50   ,  'p50b',50   ,  'p20',20    ,  'p45',45   );
        exrate      = struct( 'p_04',.04 ,  'p_01',.01  ,  'p_02',.02  ,  'p_12',.12 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utExWald(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,Case,waldmu,waldsigma,waldbarrier,exrate)
            % Computations specific to the ExWald distribution.
            testCase.ThisCase = Case;
            testCase.Dist = ExWald(waldmu,waldsigma,waldbarrier,exrate);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.DefaultParmCodes = 'rffr';
            testCase.Dist.UseSplinePDFOn(200);
            testCase.SkipMomentEst = true;  % Moments do not provide enough info to constrain parameters
            
            SetupXs(testCase,40,1000);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            testCase.CenMomentAbsTol(2) = 0.03;
            % Parameter estimation is not very good
            testCase.ParmEstAbsTol = 0.06;
            testCase.ParmEstRelTol = 0.09;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function CheckPDF(testCase)
            % Check two values supplied by Wolf Schwarz
            if testCase.ThisCase == 1
                testCase.verifyEqual(testCase.Dist.PDF(250),0.00525514318,'AbsTol',1e-6,'Does not match expected PDF value at X=250');
            elseif testCase.ThisCase == 2
                testCase.verifyEqual(testCase.Dist.PDF(250),0.003436525049,'AbsTol',1e-6,'Does not match expected PDF value at X=250');
            end
        end
        
    end
    
end  % utExWald


