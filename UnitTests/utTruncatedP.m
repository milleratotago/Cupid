classdef utTruncatedP < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
       ThisCase 
    end

    methods
        
        function testCase=utTruncatedP(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the TruncatedP distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
              case 1
                testCase.Dist = TruncatedP(Uniform(0,1),.2,.7);
                testCase.SkipEstAll = true;  % The uniform bounds don't matter as long as they are outside the cutoffs.
                testCase.Expected.Mean = (testCase.Dist.LowerBound + testCase.Dist.UpperBound)/2;
              case 2
                testCase.Dist = TruncatedP(Normal(0,1),0.11,0.93);
                testCase.EstParmCodes = 'rrff';
              case 3
                testCase.Dist = TruncatedP(Exponential(0.1),.1,0.99999);
                testCase.SkipEstML = true;  % The exponential produces a better ML at a different rate.
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            testCase.Dist.SearchOptions.MaxIter = 20000;

            SetupXs(testCase,40,5000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
            if parmCase==2
                testCase.ParmEstAbsTol(1) = .01;
                testCase.ParmEstAbsTol(2) = .01;
            end
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
                tempU = Uniform(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
                tempU = testCase.Dist.BasisRV;
                tempPDF = tempU.PDF(testCase.xvalues-testCase.Dist.LowerBound);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',2e-5,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utTruncatedP


