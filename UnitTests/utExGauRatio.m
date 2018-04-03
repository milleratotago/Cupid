classdef utExGauRatio < utContinuous;
    
    properties (ClassSetupParameter)
        parm1mu    = struct( 'p100',100 , 'p200',200 , 'p50',50   );  % normal mu
        parm2sigma = struct( 'p20',20   , 'p21',21   , 'p80',80   );  % normal sigma
        parm3ratio = struct( 'p10',10   , 'p6',6     , 'p2_8',2.8 );  % estimation is no good if parm3ratio is too small, because there is low skew & a tradeoff with normalmu
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExGauRatio(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1mu,parm2sigma,parm3ratio)
            % Computations specific to the ExGauRatio distribution.
            testCase.Dist = ExGauRatio(parm1mu,parm2sigma,parm3ratio);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
           
            testCase.MLParmTolSE = 2.5;   % ML parameter estimation is very bad
            testCase.MGFMom2RelTol = 0.02;
            testCase.ParmEstRelTol = 0.06;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExGauRatio


