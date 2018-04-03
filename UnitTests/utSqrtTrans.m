classdef utSqrtTrans < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utSqrtTrans(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the SqrtTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = SqrtTrans(Uniform(1,2));
                    % testCase.SkipAllEst = true;
                case 2
                    testCase.Dist = SqrtTrans(RNGamma(5,1/10));
                case 3
                    testCase.Dist = SqrtTrans(Exponential(0.1));
                    % testCase.SkipMLEst = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.003);
            if testCase.ThisCase==1
%                testCase.ParmEstAbsTol = .01;
            end
            if testCase.ThisCase==2
%                testCase.RawMomentAbsTol(4)=.005;
            end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
            elseif testCase.ThisCase == 2
            elseif testCase.ThisCase == 3
%                tempU = testCase.Dist.BasisRV;
%                tempPDF = tempU.PDF(testCase.xvalues-testCase.Dist.LowerBound);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utSqrtTrans


