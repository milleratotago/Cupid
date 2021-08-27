classdef utWald2MSD < utContinuous
    
    properties (ClassSetupParameter)
        Wald2MSDmu      = struct( 'p20',20 ,  'p10',10 ,  'p40',40  ,  'p50',50 );
        Wald2MSDsigma = struct( 'p50',50   ,  'p15',15 ,  'p20',20  ,  'p45',45 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utWald2MSD(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,Wald2MSDmu,Wald2MSDsigma)
            % Computations specific to the Wald2MSD distribution.
            testCase.Dist = Wald2MSD(Wald2MSDmu,Wald2MSDsigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
%            testCase.Dist.SearchOptions.MaxFunEvals = 20000;
%            testCase.Dist.SearchOptions.MaxIter = 20000;
            
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
    
end  % utWald2MSD


