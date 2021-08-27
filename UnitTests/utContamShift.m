classdef utContamShift < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utContamShift(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the ContamShift distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ContamShift(Normal(1,2),0.2,Normal(5,1.5));
                    testCase.EstParmCodes = 'fffrf';
                case 2
                    testCase.Dist = ContamShift(Beta(18,5),0.1,Normal(0.6,0.4));
                    testCase.EstParmCodes = 'ffffr';
                case 3
                    testCase.Dist = ContamShift(RNGammaMS(100,50),0.15,RNGammaMS(25,15));
                    testCase.EstParmCodes = 'frrff';
                case 4
                    testCase.Dist = ContamShift(Normal(0,1),0.1,Normal(0,0.00001));  % Constructed to give a known result.
                    testCase.EstParmCodes = 'frrff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % testCase.Dist.UseSplineCDFOn(200);
            % testCase.Dist.UseSplinePDFOn(200);
            
            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.ThisCase==2
                testCase.ParmEstAbsTol(:)=.025;
            end
            if testCase.ThisCase==4
                testCase.ParmEstRelTol(:)=.1;
                testCase.ParmEstAbsTol(:)=.02;
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
            elseif testCase.ThisCase == 4
                tempU = Normal(0,1);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utContamShift


