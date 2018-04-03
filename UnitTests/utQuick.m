classdef utQuick < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmScale = struct( 'p1',1 , 'p2',2 , 'p7',7 ,  'p8', 8  );
        parmShape = struct( 'p5',5 , 'p2',2 , 'p1',1 , 'p12',12  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utQuick(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmScale,parmShape)
            % Computations specific to the Quick distribution.
            testCase.Dist = Quick(parmScale,parmShape);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.0002);
            testCase.MGFMom2RelTol = .001;
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utQuick


