classdef utExpSum < utContinuous
    
    properties (ClassSetupParameter)
        % Remember: the rates must be different.  Testing works better if rate1 is smaller
        % because (a) the order of rates is arbitrary, and (b) PerturbParms reduces the smaller one
        % and increases the larger one.
        parmrate1  = struct( 'p_005',.005 , 'p_012',.012 , 'p_01',.01    , 'p1',1 , 'p5',5 );
        parmrate2  = struct( 'p_015',.015 , 'p_065',.065 , 'p_1',.1      , 'p2',2 , 'p12',12 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExpSum(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmrate1,parmrate2)
            % Computations specific to the ExpSum distribution.
            testCase.Dist = ExpSum(parmrate1,parmrate2);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);  % Not very accurate1
            testCase.ParmEstAbsTol = 0.01;
            testCase.ParmEstRelTol = 0.01;
            testCase.MLParmTolSE = 0.5;
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExpSum


