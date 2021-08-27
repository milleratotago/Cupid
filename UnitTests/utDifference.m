classdef utDifference < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utDifference(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Difference distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = Difference(Triangular(1,2),Uniform(0.5,1.5));
                    testCase.EstParmCodes = 'fffr';
                case 2
                    testCase.Dist = Difference(Beta(18,5),Triangular(-1,0));
                    testCase.EstParmCodes = 'fffr';
                case 3
                    testCase.Dist = Difference(RNGamma(10,2.1),Uniform(1,2));
                    testCase.EstParmCodes = 'frrf';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.Dist.UseSplinePDFOn(200);
            testCase.Dist.UseSplineCDFOn(200);
            
            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.02);
            if testCase.ThisCase==1
                testCase.RawMomentAbsTol(:)=.04;
            end
%             if testCase.ThisCase==3
testCase.ParmEstRelTol(:)=.1;
testCase.ParmEstAbsTol(:)=.1;
%                 testCase.ParmEstAbsTol(:)=.1;
%             end
            
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
        
    end % methods (Test)
    
end  % utDifference


