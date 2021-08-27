classdef utSkewNor < utContinuous
    
    
    properties (ClassSetupParameter)
        parmloc   = struct( 'n20',-20 ,  'p0',0 ,  'p40',40  ,  'p500',500 );
        parmscale = struct( 'p3',3     , 'p1',1 ,  'p50',50  ,  'p20',20    );
        parmshape = struct( 'n5',-5   ,  'n_5',-.5   ,  'p_2',.2    ,  'p4',4   );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utSkewNor(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmloc,parmscale,parmshape)
            % Computations specific to the SkewNor distribution.
            testCase.Dist = SkewNor(parmloc,parmscale,parmshape);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            
            testCase.EstParmCodes = 'frr';
            SetupXs(testCase,40,500);
            SetTolerances(testCase,0.005);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utSkewNor


