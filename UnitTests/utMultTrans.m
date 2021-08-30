classdef utMultTrans < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
       ThisCase 
    end

    methods
        
        function testCase=utMultTrans(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the MultTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
              case 1
                testCase.Dist = MultTrans(Uniform(0,1),10);
%               testCase.SkipEstAll = true;
                testCase.EstParmCodes = 'rrf';  % Don't adjust both multiplier and bounds since they can trade off.
                testCase.Expected.Mean = (testCase.Dist.LowerBound + testCase.Dist.UpperBound)/2;
              case 2
                testCase.Dist = MultTrans(Normal(100,10),0.2);
                testCase.EstParmCodes = 'rrf';  % Don't adjust both multiplier and SD since they can trade off.
              case 3
                testCase.Dist = MultTrans(Exponential(0.1),-10);
                testCase.EstParmCodes = 'rf';  % Don't adjust both rate and multiplier since they can trade off.
                testCase.SkipEstML = true;  % The exponential produces a better ML at a different rate.
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            testCase.Dist.SearchOptions.MaxIter = 20000;

            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.004);
            if testCase.ThisCase == 2
                testCase.RawMomentAbsTol(4) = 0.01;
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
            elseif testCase.ThisCase == 2
                tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
                tempU = Exponential(testCase.Dist.BasisRV.rate/testCase.Dist.Multiplier);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utMultTrans


