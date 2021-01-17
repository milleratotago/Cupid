classdef utJohnsonSU < utContinuous;
    
    properties (ClassSetupParameter)
        parm1Loc    = struct( 'p2',2  , 'p3',3    , 'n10',-10  ,  'p5',5 );
        parm2Scale  = struct( 'p2',2  , 'p4',40   , 'p3',300   ,  'p6',26 );
        parm3Alpha1 = struct( 'n5',-5 , 'n10',-10 , 'p_8',.8   ,  'p10',10   );
        parm4Alpha2 = struct( 'p2',2  , 'p3',3    , 'p2_4',2.4 ,  'p5',5 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utJohnsonSU(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1Loc,parm2Scale,parm3Alpha1,parm4Alpha2)
            % Computations specific to the JohnsonSU distribution.
            testCase.Dist = JohnsonSU(parm1Loc,parm2Scale,parm3Alpha1,parm4Alpha2);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.SearchOptions.MaxFunEvals = 20000;
            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,200);
            
            testCase.SkipEstAll = true;  % Estimation is unsatisfactory.  Estimates parameters are often quite inaccurate.
            % Set up some X values for which MLE should return approximately the true parameters:
            % npoints = 5000;
            % testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            
            testCase.ProbitNBins = 100;  % Probit estimation needs a lot of bins to recover parameters accurately.
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utJohnsonSU


