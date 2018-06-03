classdef utAttainP < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = { 1 2 3 4 };
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utAttainP(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the AttainP distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = AttainP(Normal(0.5,1),Normal(0,1),1);
                    testCase.EstParmCodes = 'rrff';
                case 2
                    testCase.Dist = AttainP(tNoncentral(20,.3),t(20),1);
                    testCase.EstParmCodes = 'frf';
                case 3
                    testCase.Dist = AttainP(FNoncentral(3,20,.1),F(3,20));
                    testCase.EstParmCodes = 'ffrff';
                case 4
                    testCase.Dist = AttainP(rNoncentral(120,-.1),r(120),-1);
                    testCase.EstParmCodes = 'frf';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            if parmCase < 4
                % Approximation is not good enough for r
                testCase.Dist.UseSplineCDFOn(200);
                testCase.Dist.UseSplinePDFOn(200);
            end

            SetupXs(testCase,40,200);
            
            testCase.Dist.SearchOptions.TolX = 1e-8;
            testCase.Dist.SearchOptions.TolFun = 1e-8;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            if testCase.ThisCase==1
                testCase.ParmEstRelTol(:)=.05;
            end
            if testCase.ThisCase==3
                testCase.ParmEstRelTol(:)=.1;
                testCase.ParmEstAbsTol(:)=.1;
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
    
end  % utAttainP


