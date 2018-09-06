classdef utBinomialZ < utDiscrete;
    
    properties (ClassSetupParameter)
        parmN   = struct('p2',2  , 'p4',4     , 'p20',20  , 'p100',100 ); % , 'p120',120, 'p130',130, 'p500',500 );
        parmp   = struct('p_1',.1, 'p_83',.83 , 'p_8',.8  , 'p_55',.55 ); % , 'p_03',.03, 'p_99',.99, 'p_73',.73 );
        parmAdj = struct('p_5',.5, 'p_25',.25 , 'p_10',.10, 'p_5b',.5 ); % , 'p_03',.03, 'p_99',.99, 'p_73',.73 );
    end
    % Note relaxed moment tolerances when using Poisson approximations
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utBinomialZ(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmp,parmAdj)

            % Computations specific to the distribution.
            testCase.Dist = BinomialZ(parmN,parmp,parmAdj);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.DiscreteX;
            
            % Set up some X values for which MLE should return close to the true parameters:
            testCase.xMLE = testCase.Dist.Random(100000,1);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            testCase.MGFh = 1.0E-4;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utBinomialZ


