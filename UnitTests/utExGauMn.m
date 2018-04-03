classdef utExGauMn < utContinuous;
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'p100',100 , 'p200',200 , 'p50',50   );  % normal mu
        parmsigma = struct( 'p20',20   , 'p21',21   , 'p80',80     );  % normal sigma
        parmmean  = struct( 'p205',205 , 'p97',97   , 'p41_3',41.3 );  % estimation is no good if parmmean is too small, because there is low skew & a tradeoff with normalmu
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExGauMn(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmsigma,parmmean)
            % Computations specific to the ExGauMn distribution.
            testCase.Dist = ExGauMn(parmmu,parmsigma,parmmean);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
           
            testCase.MLParmTolSE = 2.5;   % ML parameter estimation is very bad
            testCase.MGFMom2RelTol = 0.02;
            testCase.ParmEstRelTol = 0.06;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExGauMn


