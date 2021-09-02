classdef utConvolveFFTcIID < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = { 1 2 3 4 };
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utConvolveFFTcIID(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)

            % Adjust tolerances separately for each combination.
            
            % Computations specific to the ConvolveFFTcIID distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ConvolveFFTcIID(3,Normal(2,1));
                    testCase.EstParmCodes = 'rf';
                    SetTolerances(testCase,0.01);
                    testCase.ParmEstRelTol = 0.03;  % Moment-based estimation is not great.
                case 2
                    % The distribution of the sum of 2 Uniform(0,1) is Triangular(0,2)
                    testCase.Dist = ConvolveFFTcIID(2,Uniform(0,1));
                    testCase.EstParmCodes = 'fr';
                    SetTolerances(testCase,0.01);
                    testCase.RawMomentAbsTol(4)=.1;
                    testCase.ParmEstRelTol = 0.06;  % Moment-based estimation is not great.
                case 3
                    CommonRate = .01;
                    testCase.Dist = ConvolveFFTcIID(5,Exponential(CommonRate));
                    testCase.EstParmCodes = 'r';
                    SetTolerances(testCase,0.01);
                    testCase.RawMomentRelTol(1) = 0.01;
                    testCase.RawMomentAbsTol(4)=.1;
                    testCase.CenMomentAbsTol(2) = 0.5;  % relative to mean of 500
                    testCase.ParmEstRelTol = 0.02;
                case 4
                    testCase.Dist = ConvolveFFTcIID(3,RNGamma(2,.01));
                    testCase.EstParmCodes = 'fr';
                    SetTolerances(testCase,0.01);
                    testCase.ParmEstRelTol = 0.02;
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
                tempU = RNGamma(testCase.Dist.NDists,testCase.Dist.BasisRV.rate);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 4
                % Sum of identical RNGammas is RNGamma with sum of k's
                tempU = RNGamma(testCase.Dist.NDists*testCase.Dist.BasisRV.N,testCase.Dist.BasisRV.Rate);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utConvolveFFTcIID


