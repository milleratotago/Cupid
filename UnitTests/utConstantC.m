classdef utConstantC < utContinuous;
    
    properties (ClassSetupParameter)
        thisconstant    = struct( 'p10',10 , 'p50',50 , 'p200',200 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utConstantC(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,thisconstant)
            % Computations specific to the distribution.
            testCase.Dist = ConstantC(thisconstant);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,2000);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = thisconstant;
            testCase.Expected.SD = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

%            testCase.SkipMGFs = true;
            testCase.SkipAllEst = true;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utConstantC


