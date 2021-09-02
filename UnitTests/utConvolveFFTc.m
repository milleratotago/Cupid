classdef utConvolveFFTc < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = { 1 2 3 4 5 };
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utConvolveFFTc(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)

            % Adjust tolerances separately for each combination.
            
            % Computations specific to the ConvolveFFTc distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ConvolveFFTc(Normal(0.5,1),Normal(0,1));
                    testCase.EstParmCodes = 'rfff';
                    SetTolerances(testCase,0.01);
                case 2
                    % The distribution of the sum of 2 Uniform(0,1) is Triangular(0,2)
                    testCase.Dist = ConvolveFFTc(Uniform(0,1),Uniform(0,1));
                    testCase.EstParmCodes = 'rrff';
                    SetTolerances(testCase,0.01);
                    testCase.RawMomentAbsTol(4)=.1;
                    testCase.ParmEstRelTol = 0.02;
                case 3
                    CommonRate = .01;
                    testCase.Dist = ConvolveFFTc(RNGamma(4,CommonRate),Exponential(CommonRate));
                    testCase.EstParmCodes = 'rff';
                    SetTolerances(testCase,0.01);
                    testCase.RawMomentRelTol(1) = 0.01;
                    testCase.CenMomentRelTol(1) = 0.01;
                    testCase.CenMomentAbsTol(2) = 0.1;
                    testCase.ParmEstRelTol = 0.02;
                case 4
                    testCase.Dist = ConvolveFFTc(Normal(200,20),Exponential(.01));
                    testCase.EstParmCodes = 'frf';
                    SetTolerances(testCase,0.01);
                    testCase.ParmEstRelTol = 0.02;
                    testCase.CenMomentAbsTol(2) = 0.03;
                case 5  % same as case 4 but with 3 distributions
                    testCase.Dist = ConvolveFFTc(Normal(100,20/sqrt(2)),Normal(100,20/sqrt(2)),Exponential(.01));
                    testCase.EstParmCodes = 'frf';
                    SetTolerances(testCase,0.01);
                    testCase.ParmEstRelTol = 0.02;
                    testCase.CenMomentAbsTol(2) = 0.03;
                    testCase.SkipEstAll = true;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

%            testCase.Dist.UseSplineCDFOn(200);
            % testCase.Dist.UseSplinePDFOn(1000);  % ConvolveFFTc uses its own spline approximation.
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
                 tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
                 tempPDF = tempU.PDF(testCase.xvalues);
                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-4,'RelTol',1e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
                % Sum of two uniforms is triangular
                tempU = Triangular(testCase.Dist.Minimum,testCase.Dist.Maximum);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',2e-4,'RelTol',2e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
                % Sum of Gamma and same-rate exponential is another gamma with increased N
                tempU = RNGamma(testCase.Dist.BasisRV{1}.N+1,testCase.Dist.BasisRV{1}.Rate);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 4
                % Sum of normal and exponential is ExGaussian
                tempU = ExGauss(testCase.Dist.BasisRV{1}.Mean,testCase.Dist.BasisRV{1}.SD,testCase.Dist.BasisRV{2}.rate);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 5
                % Sum of normal and exponential is ExGaussian
                NorMn = testCase.Dist.BasisRV{1}.Mean + testCase.Dist.BasisRV{2}.Mean;
                NorSD = sqrt(testCase.Dist.BasisRV{1}.Variance + testCase.Dist.BasisRV{2}.Variance);
                tempU = ExGauss(NorMn,NorSD,testCase.Dist.BasisRV{3}.rate);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utConvolveFFTc


