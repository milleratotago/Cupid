classdef utExWaldMSM < utContinuous;
    
    
    properties (ClassSetupParameter)
        Case       = struct( 'case1',1  , 'case2',2  , 'case3',3  , 'case4',4  );
        p1WaldMean = struct( 'p200',200 , 'p300',300 , 'p500',500 , 'p750',750 );
        p2WaldSD   = struct( 'p50',50   , 'p75',75   , 'p90',90   , 'p145',145 );
        p3ExMean   = struct( 'p205',205 , 'p100',100 , 'p150',50  , 'p200',200 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utExWaldMSM(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,Case,p1WaldMean,p2WaldSD,p3ExMean)
            % Computations specific to the ExWaldMSM distribution.
            testCase.ThisCase = Case;
            testCase.Dist = ExWaldMSM(p1WaldMean,p2WaldSD,p3ExMean);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.DefaultParmCodes = 'rfr';
            testCase.Dist.UseSplinePDFOn(200);
            testCase.SkipMomentEst = true;  % Moments do not provide enough info to constrain parameters
            
            SetupXs(testCase,40,1000);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            % Parameter estimation is not very good
            testCase.ParmEstAbsTol = 0.06;
            testCase.ParmEstRelTol = 0.09;
            testCase.RawMomentRelTol(4) = .03;            
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utExWaldMSM


