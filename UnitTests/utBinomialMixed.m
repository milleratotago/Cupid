classdef utBinomialMixed < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = { 1 2 3 4 };
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end

    methods
        
        function testCase=utBinomialMixed(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    ThesePs = 0.4*ones(10,1);
                case 2
                    ThesePs = 0.6*ones(25,1);
                case 3
                    ThesePs = (1:10)/20;
                case 4
                    ThesePs = (2.5:5:97.5)/100;
            end

            testCase.Dist = BinomialMixed(ThesePs);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.LowerBound:testCase.Dist.UpperBound;
            
            % Set up some X values for which MLE should return the true parameters:
            npoints = 500;
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            testCase.SkipEstAll = true;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
                 tempU = Binomial(10,.4);
                 tempPDF = tempU.PDF(testCase.xvalues);
                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-4,'RelTol',1e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
                tempU = Binomial(25,.6);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',2e-4,'RelTol',2e-4,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utBinomialMixed


