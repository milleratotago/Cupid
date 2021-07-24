classdef utContamStretch < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utContamStretch(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the ContamStretch distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    StretchDist = Triangular(1.3,1.5);
                    testCase.EstParmCodes = 'ffrff';
                    testCase.Dist = ContamStretch(RNGamma(4,0.2),0.15,StretchDist);
                case 2
                    StretchDist = Triangular(1.1,1.35);
                    testCase.Dist = ContamStretch(Beta(18,5),0.1,StretchDist);
                    testCase.EstParmCodes = 'frfff';
                case 3
                    StretchDist = Uniform(1.1,1.35);
                    testCase.Dist = ContamStretch(RNGammaMS(100,50),0.15,StretchDist);
                    testCase.EstParmCodes = 'frrff';
                case 4
                    StretchDist = LognormalMS(1.8,0.12);
                    testCase.Dist = ContamStretch(Weibull(10,1.2,13),0.05,StretchDist);
                    testCase.EstParmCodes = 'fffrff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % testCase.Dist.UseSplineCDFOn(200);
            % testCase.Dist.UseSplinePDFOn(200);
            
            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.ThisCase==4
                SetTolerances(testCase,0.01);
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
%                 tempU = Normal(0,1);
%                 tempPDF = tempU.PDF(testCase.xvalues-100);
%                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utContamStretch


