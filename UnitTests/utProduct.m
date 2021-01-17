classdef utProduct < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utProduct(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Product distribution.
            testCase.ThisCase  = parmCase;
            testCase.SkipEstAll = true;  % ESTIMATION IS TOO SLOW
            switch parmCase
                case 1
                    testCase.Dist = Product(RNGamma(10,.1),RNGamma(10,10)); % Lognormal(1,.1),Uniform(1,2));
                    % testCase.EstParmCodes = 'rfff';
                case 2
                    testCase.Dist = Product(Beta(5,18),Wald2(.25,20));
                    % testCase.EstParmCodes = 'ffrf';
                case 3
                    testCase.Dist = Product(RNGamma(10,.1),Uniform(1,2));
                    % testCase.EstParmCodes = 'rfff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.UseSplineCDFOn(1000);  % Very, very slow without splines.
            testCase.Dist.UseSplinePDFOn(1000);
            
            SetupXs(testCase,80,500);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            testCase.ParmEstRelTol(:) = .1;
            testCase.ParmEstAbsTol(:) = .1;
            
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
    
end  % utProduct


