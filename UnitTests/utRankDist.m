classdef utRankDist < utDiscrete;
    
    properties (ClassSetupParameter)
        parmCase = {1 2 3 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end

    methods
        
        function testCase=utRankDist(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the RankDist distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = RankDist(10,Uniform(0,1),Normal(0,1));
                    testCase.EstParmCodes = 'fffrr';
                case 2
                    testCase.Dist = RankDist(5,Exponential(1),Triangular(0,5));
                    testCase.EstParmCodes = 'frff';
                case 3
                    testCase.Dist = RankDist(4,Poisson(20),Normal(20,2));
                    testCase.EstParmCodes = 'frff';
                case 4
                    testCase.Dist = RankDist(8,Normal(20,5),Poisson(16));
                    testCase.EstParmCodes = 'frff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.XsToPlot;
            testCase.xMLE = testCase.Dist.Random(100000,1);
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utRankDist


