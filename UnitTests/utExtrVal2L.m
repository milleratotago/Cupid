classdef utExtrVal2L < utContinuous;
    
    properties (ClassSetupParameter)
        parmscale = struct( 'p5',5   , 'p300',300, 'p22',22   , 'p84',84 , 'p2000',2000 );
        parmshape = struct( 'p5',5   , 'p10',10  , 'p3',3     , 'p4',4   , 'p11',11 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExtrVal2L(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmscale,parmshape)
            % Computations specific to the ExtrVal2L distribution.
            testCase.Dist = ExtrVal2L(parmscale,parmshape);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
%            testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExtrVal2L


