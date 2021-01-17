classdef utOrderDisc < utDiscrete
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utOrderDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Order distribution.
            testCase.ThisCase  = parmCase;
%             testCase.SkipEstAll = true;  % Parameter estimation is slow
            switch parmCase
                case 1
                    testCase.Dist = Order(3,UniformInt(1,10),Poisson(10),Binomial(20,.5));
                    testCase.EstParmCodes = 'fffffr';
                case 2
                    testCase.Dist = Order(2,Geometric(.1),Poisson(10));
                    testCase.EstParmCodes = 'frf';
                case 3
                    testCase.Dist = Order(3,Binomial(10,0.1),Binomial(20,0.1),Binomial(30,0.1));
                    testCase.EstParmCodes = 'ffrffff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

%            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
%            testCase.Dist.SearchOptions.MaxIter = 20000;
            
            SetupXs(testCase);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
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
%                tempU = Binomial(testCase.Dist.BasisRV{1}.rate*3);
%                tempPDF = tempU.PDF(testCase.xvalues);
%                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utOrderDisc


