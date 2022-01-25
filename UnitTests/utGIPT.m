classdef utGIPT < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1, 2, 3, 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utGIPT(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the GIPT distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = GIPT(Normal(0.5,1),Uniform(eps,1-eps));
                    testCase.EstParmCodes = 'rrff';
                case 2
                    testCase.Dist = GIPT(Normal(0.5,1),Beta(10,10));
                    testCase.EstParmCodes = 'rrff';
                case 3
                    testCase.Dist = GIPT(Exponential(1/300),Beta(30,10));
                    testCase.EstParmCodes = 'rff';
                case 4
                    testCase.Dist = GIPT(RNGammaMS(750,126),Beta(10,30));
                    testCase.EstParmCodes = 'rrff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            testCase.Dist.SearchOptions.TolX = 1e-8;
            testCase.Dist.SearchOptions.TolFun = 1e-8;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
                                tempN = Normal(testCase.Dist.X.Mean,testCase.Dist.X.SD);
                                tempPDF = tempN.PDF(testCase.xvalues);
                                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
            elseif testCase.ThisCase == 3
            end
        end
        
    end
    
end  % utGIPT


