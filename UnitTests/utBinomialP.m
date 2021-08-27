classdef utBinomialP < utDiscrete
    
    properties (ClassSetupParameter)
        parmN = struct('p2',2  , 'p4',4     , 'p20',20, 'p100',100 , 'p120',120, 'p130',130, 'p500',500 );
        parmp = struct('p_1',.1, 'p_83',.83 , 'p_8',.8, 'p_55',.55 , 'p_03',.03, 'p_99',.99, 'p_73',.73 );
    end
    % Note relaxed moment tolerances when using Poisson approximations
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utBinomialP(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmp)
            % Computations specific to the distribution.
            testCase.Dist = BinomialP(parmN,parmp);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.LowerBound:(1/parmN):testCase.Dist.UpperBound;
            
            % Set up some X values for which MLE should return the true parameters:
            npoints = 500;
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmp;
            testCase.Expected.Variance = parmp*(1-parmp)/parmN;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.Dist.Approx>=2  % Poisson approximation moments are not very accurate
                testCase.RawMomentRelTol(3) = 0.04;
                testCase.CenMomentRelTol(3) = 0.04;
                testCase.KurtRelTol = 0.04;
            end

%            if parmN==100   % MGF not very accurate in these cases.
%                testCase.MGFMom2AbsTol = .01;
%            elseif parmN==500
%                testCase.MGFMom2RelTol = .10;
%            end

            testCase.MGFh = 1.0E-4;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utBinomialP


