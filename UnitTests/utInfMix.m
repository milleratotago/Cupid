classdef utInfMix < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utInfMix(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the InfMix distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = InfMix(Normal(0,1),Uniform(0.9999,1.0001),2);  % Constructed to give a known result.
                    testCase.Dist.DefaultParmCodes = 'rfff';
                    testCase.SkipAllEst = true;  % Skip in interests of speed
                case 2
                    testCase.Dist = InfMix(Normal(.09091,1),Beta(5,50),1);
                    testCase.Dist.DefaultParmCodes = 'frff';
                    testCase.SkipAllEst = true;  % Skip in interests of speed
                case 3
                    testCase.Dist = InfMix(Exponential(1),Triangular(0.5,1.5),1);
                    testCase.Dist.DefaultParmCodes = 'frf';
                    testCase.SkipAllEst = true;  % Skip in interests of speed
                case 4
                    testCase.Dist = InfMix(RNGamma(5,.1),Cosine(5,3),1);
                    testCase.Dist.DefaultParmCodes = 'frff';  % DO THIS ESTIMATION JUST AS ONE EXAMPLE.
                case 5
                    testCase.Dist = InfMix(RNGamma(5,.1),Cosine(5,2),1);
                    testCase.Dist.DefaultParmCodes = 'ffrf';  % Note adjusting parameter of ParmDist rather than ParentDist
                    testCase.SkipAllEst = true;  % Skip in interests of speed
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % testCase.Dist.UseSplineCDFOn(200);
            % testCase.Dist.UseSplinePDFOn(200);
            
            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
%             if testCase.ThisCase==1
%                 testCase.ParmEstRelTol(:)=.05;
%             end
%             if testCase.ThisCase==3
%                 testCase.ParmEstRelTol(:)=.1;
%                 testCase.ParmEstAbsTol(:)=.1;
%             end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
                tempU = Normal(0,1);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
            elseif testCase.ThisCase == 3
            elseif testCase.ThisCase == 4
            end
        end
        
    end
    
end  % utInfMix


