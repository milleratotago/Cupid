classdef utExGauss < utContinuous
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'p100',100   , 'p200',200 , 'p500',500   );  % normal mu
        parmsigma = struct( 'p20',20     , 'p21',21   , 'p80',80     );  % normal sigma
        parmrate  = struct( 'p_005',.005 , 'p_01',.01 , 'p_025',.025 );  % exp mean = 200, 100, 40
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExGauss(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmsigma,parmrate)
            % Computations specific to the ExGauss distribution.
            testCase.Dist = ExGauss(parmmu,parmsigma,parmrate);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
           
            testCase.MLParmTolSE = 2.5;   % ML parameter estimation is very bad
            testCase.MGFMom2RelTol = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExGauss


