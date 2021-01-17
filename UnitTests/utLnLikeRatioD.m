classdef utLnLikeRatioD < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utLnLikeRatioD(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LikeRatio distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LnLikeRatioD(Binomial(10,.5),Normal(5,1),Normal(5.5,1));
                    testCase.EstParmCodes = 'frffff';
%                     testCase.SkipEstAll = true;  % ESTIMATION IS SLOW SO ONLY DO ONE
                case 2
                    testCase.Dist = LnLikeRatioD(NegativeBinomial(5,.5),Normal(5,5),Normal(5.5,5));
                    testCase.EstParmCodes = 'frffff';
                    testCase.Dist.XGrain = 1000;  % Numerical problems with X values
%                     testCase.SkipEstAll = true;  % ESTIMATION IS SLOW SO ONLY DO ONE
                case 3
                    testCase.Dist = LnLikeRatioD(Poisson(10),Normal(8,3.2),Normal(10,3.2));
                    testCase.EstParmCodes = 'rffff';
%                     testCase.SkipEstAll = true;  % ESTIMATION IS SLOW SO ONLY DO ONE
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            %             testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            %             testCase.Dist.SearchOptions.MaxIter = 20000;
            
%             testCase.Dist.UseSplineCDFOn(200);  % Already using splines

            testCase.SetupXs;

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.ThisCase==1
                testCase.ParmEstRelTol(:)=.05;
            end
            if testCase.ThisCase==2
                testCase.ParmEstRelTol(:)=.2;
                testCase.ParmEstAbsTol(:)=.2;
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
    
end  % utLnLikeRatioD



