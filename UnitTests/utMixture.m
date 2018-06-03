classdef utMixture < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utMixture(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Mixture distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = Mixture(0.2,Triangular(1,2),Uniform(0.5,1.5));
                    testCase.EstParmCodes = 'ffffr';
                case 2
                    testCase.Dist = Mixture(0.5,Beta(18,5),Beta(21,2));
                    testCase.EstParmCodes = 'ffffr';
                case 3
                    testCase.Dist = Mixture(0.8,RNGamma(10,.1),Triangular(50,150));
                    testCase.EstParmCodes = 'frrff';
                case 4
                    testCase.Dist = Mixture(0.5,Normal(0,1),Normal(0,1));  % Constructed to give a known result.
                    testCase.EstParmCodes = 'frrff';
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
            elseif testCase.ThisCase == 2
            elseif testCase.ThisCase == 3
            elseif testCase.ThisCase == 4
                tempU = Normal(0,1);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utMixture


