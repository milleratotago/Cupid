classdef utMonotoneTransDisc < utDiscrete
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utMonotoneTransDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the MonotoneTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    f = @(x) 2*x+2;
                    finv = @(x) (x-2)/2;
                    testCase.Dist = MonotoneTrans(UniformInt(0,10),f,finv);
                    testCase.Dist.XGrain = 1000;
                    testCase.SkipEstMoment = true;    % Not good at estimating these int parms
                    testCase.SkipEstProbitAll = true;
                case 2
                    f = @(x) x/2+3;
                    finv = @(x) (x-3)*2;
                    testCase.Dist = MonotoneTrans(Binomial(88,.33),f,finv);
                    testCase.EstParmCodes = 'fr';
                    testCase.Dist.XGrain = 1000;
                case 3
                    f = @(x) 2*x+2;
                    finv = @(x) (x-2)/2;
                    testCase.Dist = MonotoneTrans(Poisson(20.1),f,finv);
                    testCase.Dist.XGrain = 1000;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
%            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
%            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.02);
%            if testCase.ThisCase==1
%                testCase.ParmEstAbsTol = .01;
%            end
%            if testCase.ThisCase==2
%                testCase.RawMomentAbsTol(4)=.005;
%            end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
%                tempU = UniformInt(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
%                tempPDF = tempU.PDF(testCase.xvalues);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
%                tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
%                tempPDF = tempU.PDF(testCase.xvalues);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
%                tempU = testCase.Dist.BasisRV;
%                tempPDF = tempU.PDF(testCase.Dist.InverseFunc(testCase.xvalues))/2;
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utMonotoneTransDisc


