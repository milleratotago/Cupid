classdef utLnLikeRatioC < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utLnLikeRatioC(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the LikeRatio distribution.
            NBins = 300;
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = LnLikeRatioC(Beta(4,4),Exponential(1),Exponential(0.5),NBins);
                    testCase.EstParmCodes = 'frff';
                case 2
                    testCase.Dist = LnLikeRatioC(Beta(4.3,2),Beta(2,2),Beta(4,2),NBins);
                    testCase.EstParmCodes = 'rrffff';
                case 3
                    testCase.Dist = LnLikeRatioC(TruncatedX(RNGamma(10,.1),1,500),RNGamma(8,.1),RNGamma(9,.1),NBins);
                    testCase.EstParmCodes = 'frffffff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            %             testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            %             testCase.Dist.SearchOptions.MaxIter = 20000;
            
%             testCase.Dist.UseSplineCDFOn(200);  % Already using splines
            testCase.Dist.UseSplinePDFOn(200);

            SetupXs(testCase,40,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.ThisCase==1
                testCase.ParmEstRelTol(:)=.05;
            end
            if testCase.ThisCase==2
                testCase.ParmEstRelTol(:)=.2;
                testCase.ParmEstAbsTol(:)=.2;
                testCase.RawMomentAbsTol(1) = 0.02;  % Bad numerical problems with this one
                testCase.CenMomentAbsTol(1) = 0.02;
                testCase.CenMomentAbsTol(2) = 0.01;
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
    
end  % utLnLikeRatioC

