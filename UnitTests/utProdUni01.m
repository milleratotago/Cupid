classdef utProdUni01 < utContinuous
    
    properties (ClassSetupParameter)
        parmK   = struct( 'K1',1 , 'K2',2  , 'K3',3    , 'K5',5  );  % Big numerical problems already at K = 10
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utProdUni01(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmK)

            testCase.Dist = ProdUni01(parmK);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
%            testCase.SkipEstMoment = true;  % Moments do not give much info about df's
%            testCase.SkipEstML = true;  % I can't find X values for which the MLE is the current dfs.
            % Set up some X values for which MLE should return (very close to) the true parameters:
            % Skipped because I can'F find X values for which the MLE is the current df.
%             npoints = 2000; % 20000;
%             testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);

            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utProdUniform1


