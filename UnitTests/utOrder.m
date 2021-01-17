classdef utOrder < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utOrder(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Order distribution.
            testCase.ThisCase  = parmCase;
%             testCase.SkipEstAll = true;  % Parameter estimation is slow
            switch parmCase
                case 1
%                     testCase.Dist = Order(3,Uniform(0,2),Normal(0,2),Uniform(1,3));  % Very slow
                    testCase.Dist = Order(3,Uniform(0,2),Triangular(0,2),TriangularG(0,.5,2));
                    testCase.EstParmCodes = 'ffffffrf';
                case 2
                    testCase.Dist = Order(2,Normal(100,10),Gamma(10,.1));
                    testCase.EstParmCodes = 'frfff';
                    testCase.Dist.UseSplinePDF = true;   % Speeds things up enormously
                case 3
                    testCase.Dist = Order(1,Exponential(0.1),Exponential(0.1),Exponential(0.1));
                    testCase.EstParmCodes = 'frff';
%                     testCase.Dist.UseSplinePDF = true;   % Speeds things up enormously
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

%            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
%            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.02);
%            if testCase.ThisCase==1
%                testCase.ParmEstAbsTol = .01;
%            end
%            if testCase.ThisCase==2
%                testCase.RawMomentAbsTol(4)=.005;
%            end
            
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
                tempU = Exponential(testCase.Dist.BasisRV{1}.rate*3);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utOrder


