classdef utLogitLogistic < utContinuous
% Special class to test combination of logistic & logit transformations,
% which should result in the original X
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
        BasisRVs
    end
    
    methods
        
        function testCase=utLogitLogistic(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
            testCase.BasisRVs{1} = Normal(0.5,0.5);
            testCase.BasisRVs{2} = Beta(22,5);
            testCase.BasisRVs{3} = RNGamma(4,10);
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LogitTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LogitTrans(LogisticTrans(testCase.BasisRVs{1}));
                    %               testCase.SkipEstAll = true;
                case 2
                    testCase.Dist = LogitTrans(LogisticTrans(testCase.BasisRVs{2}));
                case 3
                    testCase.Dist = LogitTrans(LogisticTrans(testCase.BasisRVs{3}));
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.ThisCase==2
                testCase.RawMomentAbsTol(4)=.01;
            end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            tempDist = testCase.BasisRVs{testCase.ThisCase};
            tempPDF = tempDist.PDF(testCase.xvalues);
            testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
        end
        
    end
    
end  % utLogitLogistic


