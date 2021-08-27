classdef utExpLogDisc < utDiscrete
% Special class to test combination of Log & Exp transformations,
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
        
        function testCase=utExpLogDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
            testCase.BasisRVs{1} = Binomial(9,.44);
            testCase.BasisRVs{2} = UniformInt(0,5);
            testCase.BasisRVs{3} = Poisson(6.55);
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the ExpTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ExpTrans(LogTrans(testCase.BasisRVs{1}));
                    testCase.EstParmCodes = 'fr';
                case 2
                    testCase.Dist = ExpTrans(LogTrans(testCase.BasisRVs{2}));
%                     testCase.EstParmCodes = 'r';
                    testCase.SkipEstAll = true;
                case 3
                    testCase.Dist = ExpTrans(LogTrans(testCase.BasisRVs{3}));
                    testCase.EstParmCodes = 'r';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            testCase.Dist.XGrain = 1000;  % Be very tolerant of numerical errors due to the double non-linear conversion.
            testCase.Dist.SetBinEdges;    % XGrain only has an effect when bin edges are set, so re-do.

            testCase.SetupXs;  % xvalues at which PDF, CDF, etc should be evaluated & xMLE used to estimate MLE
            
            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
%             if testCase.ThisCase==2
%                 testCase.RawMomentAbsTol(4)=.01;
%             end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            tempDist = testCase.BasisRVs{testCase.ThisCase};
            tempDist.XGrain = 1000;  % Be very tolerant of numerical errors (i.e. in the .xvalues) due to the double non-linear conversion.
            tempDist.SetBinEdges;    % XGrain only has an effect when bin edges are set, so re-do.
            tempPDF = tempDist.PDF(testCase.xvalues);
            testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
        end
        
    end
    
end  % utLogExpDisc


