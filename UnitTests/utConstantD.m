classdef utConstantD < utDiscrete;
    
    properties (ClassSetupParameter)
        thisConstant = struct('n2',-2  , 'p4',4, 'n20',-20, 'p100',100, 'p500',500 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utConstantD(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,thisConstant)
            % Computations specific to the distribution.
            testCase.Dist = ConstantD(thisConstant);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = thisConstant;
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = thisConstant;
            testCase.Expected.Variance = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

%            testCase.SkipMGFs = true;
            testCase.SkipAllEst = true;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utConstantD


