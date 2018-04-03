classdef utBeta < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmA = struct( 'p1',1 , 'p2',2 , 'p7',7 ,  'p8', 8  );
        parmB = struct( 'p5',5 , 'p2',2 , 'p1',1 , 'p12',12  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utBeta(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmA,parmB)
            % Computations specific to the Beta distribution.
            testCase.Dist = Beta(parmA,parmB);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            % testCase.Dist.SearchOptions.TolFun = 1e-14;
            % testCase.Dist.SearchOptions.TolX = 1e-14;

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.0002);
            testCase.MGFMom2AbsTol = 0.005;
            testCase.MGFMom2RelTol = 0.005;
            testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.01;
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.01;
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utBeta


