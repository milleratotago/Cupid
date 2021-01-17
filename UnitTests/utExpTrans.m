classdef utExpTrans < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utExpTrans(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the ExpTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ExpTrans(Uniform(1,2));
                    % testCase.SkipEstAll = true;
                case 2
                    testCase.Dist = ExpTrans(Normal(1,1));
                case 3
                    testCase.Dist = ExpTrans(Exponential(5));
                    % testCase.SkipEstML = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            
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
%                tempU = Uniform(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
%                tempPDF = tempU.PDF(testCase.xvalues);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
                tempU = Lognormal(testCase.Dist.BasisRV.Mean,testCase.Dist.BasisRV.SD);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
%                tempU = testCase.Dist.BasisRV;
%                tempPDF = tempU.PDF(testCase.xvalues-testCase.Dist.LowerBound);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utExpTrans


