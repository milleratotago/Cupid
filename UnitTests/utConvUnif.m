classdef utConvUnif < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = { 1 2 3 4 };
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utConvUnif(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)

            % Adjust tolerances separately for each combination.
            
            % Computations specific to the ConvUnif distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ConvUnif(Normal(0.5,1),-1,1);
                    testCase.EstParmCodes = 'rfff';
                    SetTolerances(testCase,0.01);
                case 2
                    % The distribution of the sum of 2 Uniform(0,1) is Triangular(0,2)
                    testCase.Dist = ConvUnif(Uniform(0,1),0,1);
                    testCase.EstParmCodes = 'rrff';
                    SetTolerances(testCase,0.01);
                    testCase.RawMomentAbsTol(4)=.1;
                    testCase.ParmEstRelTol = 0.02;
                case 3
                    CommonRate = .01;
                    testCase.Dist = ConvUnif(RNGamma(4,CommonRate),100,200);
                    testCase.EstParmCodes = 'rrff';
                    SetTolerances(testCase,0.01);
                    testCase.RawMomentRelTol(1) = 0.01;
                    testCase.CenMomentRelTol(1) = 0.01;
                    testCase.ParmEstRelTol = 0.02;
                case 4
                    testCase.Dist = ConvUnif(Normal(200,20),-5,7);
                    testCase.EstParmCodes = 'rrff';
                    SetTolerances(testCase,0.01);
                    testCase.ParmEstRelTol = 0.02;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

%            testCase.Dist.UseSplineCDFOn(200);
%            testCase.Dist.UseSplinePDFOn(1000);  % Use a lot of points to get a good approximation.
            testCase.Dist.SearchOptions.TolX = 1e-8;
            testCase.Dist.SearchOptions.TolFun = 1e-8;
            
            SetupXs(testCase,40,2000);
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
                 tempU = Convolution(Normal(0.5,1),Uniform(-1,1));
                 tempPDF = tempU.PDF(testCase.xvalues);
                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-4,'RelTol',1e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
                % Sum of two uniforms is triangular
                tempU = Triangular(testCase.Dist.Minimum,testCase.Dist.Maximum);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',2e-4,'RelTol',2e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
                tempU = Convolution(testCase.Dist.BasisRV1,Uniform(100,200));
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 4
                tempU = Convolution(Normal(200,20),Uniform(-5,7));
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utConvUnif


