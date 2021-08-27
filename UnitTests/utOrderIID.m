classdef utOrderIID < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utOrderIID(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the OrderIID distribution.
            testCase.ThisCase  = parmCase;
%             testCase.SkipEstAll = true;  % 6 minutes skipping these, not parallel.
            switch parmCase
                case 1
                    testCase.Dist = OrderIID(1,3,Uniform(0,2));
                case 2
                    testCase.Dist = OrderIID(2,5,Normal(100,10));
                    testCase.Dist.UseSplinePDF = true;   % Speeds things up enormously
                case 3
                    testCase.Dist = OrderIID(1,10,Exponential(0.1));
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

%            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
%            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,1000);
            
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
%                tempU = Uniform(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
%                tempPDF = tempU.PDF(testCase.xvalues);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
%                tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
%                tempPDF = tempU.PDF(testCase.xvalues);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
                tempU = Exponential(testCase.Dist.BasisRV.rate*testCase.Dist.SampleSize);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utOrderIID


