classdef utMinBound < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utMinBound(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the MinBound distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = MinBound(Uniform(10,20),Uniform(10,20));
                    testCase.Dist.DefaultParmCodes = 'ffrr';
                case 2
                    testCase.Dist = MinBound(Beta(18,5),Beta(21,2));
                    testCase.Dist.DefaultParmCodes = 'fffr';
                case 3
                    testCase.Dist = MinBound(RNGamma(10,10),Uniform(1,2));
                    testCase.Dist.DefaultParmCodes = 'rrff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.UseSplineCDFOn(1000);
            testCase.Dist.UseSplinePDFOn(1000);
            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,80,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.ParmEstRelTol = 10 * testCase.ParmEstRelTol;
            
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
                %                tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
                %                tempPDF = tempU.PDF(testCase.xvalues);
                %                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
                %                tempU = testCase.Dist.BasisRV;
                %                tempPDF = tempU.PDF(testCase.xvalues-testCase.Dist.LowerBound);
                %                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utMinBound


