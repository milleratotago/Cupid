classdef utHypergeometric < utDiscrete;
    
    properties (ClassSetupParameter)
        parmN = struct('N28',28, 'N30',30, 'N49',49, 'N51',51, 'N56',56);
        parmK = struct('K14',14, 'K16',16, 'K25',25, 'K36',36, 'K17',17);
        parmn = struct('n14',14, 'n22',22, 'n35',35, 'n20',20, 'n30',30);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utHypergeometric(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmK,parmn)
            % Computations specific to the distribution.
            testCase.Dist = Hypergeometric(parmN,parmK,parmn);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.SetupXs;
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Expect some problems due to truncation of discrete distribution

            testCase.SkipEstAll = true;
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utHypergeometric


