classdef utLinearTrans < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utLinearTrans(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LinearTrans distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LinearTrans(Uniform(0,1),2,10);
                    testCase.Dist.DefaultParmCodes = 'rrff';  % Don't adjust both constants and bounds since they can trade off.
                    %               testCase.SkipAllEst = true;
                    testCase.Expected.Mean = (testCase.Dist.LowerBound + testCase.Dist.UpperBound)/2;
                case 2
                    testCase.Dist = LinearTrans(Normal(0,1),100,10);
                    testCase.Dist.DefaultParmCodes = 'rrff';  % Don't adjust both constants and mean/sd since they can trade off.
                case 3
                    testCase.Dist = LinearTrans(Exponential(0.1),2,5);
                    testCase.Dist.DefaultParmCodes = 'rff';  % Don't adjust both constants and rate since they can trade off.
                    testCase.SkipMLEst = true;  % The exponential produces a better ML at a different rate.
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
            if testCase.ThisCase==1
                testCase.ParmEstAbsTol = .01;
            end
            if testCase.ThisCase==2
                testCase.RawMomentAbsTol(4)=.005;
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
                tempU = testCase.Dist.BasisRV;
                tempPDF = tempU.PDF((testCase.xvalues-testCase.Dist.LowerBound)/2)/2;
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utLinearTrans


