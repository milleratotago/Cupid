classdef utConvolutionDisc < utDiscrete;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = { 1 2 3 4 };
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utConvolutionDisc(varargin)  % Constructor
            testCase=testCase@utDiscrete(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the Convolution distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = Convolution(Binomial(10,.5),Poisson(2.3));
                    testCase.EstParmCodes = 'frf';
                case 2
                    testCase.Dist = Convolution(UniformInt(0,10),Geometric(.10));
                    testCase.EstParmCodes = 'ffr';
                case 3
                    testCase.Dist = Convolution(Poisson(2.5),Poisson(3.2));
                    testCase.EstParmCodes = 'rf';  % Parameters trade off if both varied.
                case 4
                    testCase.Dist = Convolution(Binomial(10,.5),Binomial(4,.9));
                    testCase.EstParmCodes = 'frff';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            
            %testCase.Dist.SearchOptions.TolX = 1e-8;
            %testCase.Dist.SearchOptions.TolFun = 1e-8;
            
            testCase.SetupXs;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.ThisCase == 1
%                  tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
%                  tempPDF = tempU.PDF(testCase.xvalues);
%                  testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-4,'RelTol',1e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 2
%                 % Sum of two uniforms is triangular
%                 tempU = Triangular(testCase.Dist.Minimum,testCase.Dist.Maximum);
%                 tempPDF = tempU.PDF(testCase.xvalues);
%                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',2e-4,'RelTol',2e-4,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 3
%                 % Sum of Gamma and same-rate exponential is another gamma with increased N
%                 tempU = RNGamma(testCase.Dist.BasisRV1.N+1,testCase.Dist.BasisRV2.rate);
%                 tempPDF = tempU.PDF(testCase.xvalues);
%                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            elseif testCase.ThisCase == 4
%                 % Sum of normal and exponential is ExGaussian
%                 tempU = ExGauss(testCase.Dist.BasisRV1.Mean,testCase.Dist.BasisRV1.SD,testCase.Dist.BasisRV2.rate);
%                 tempPDF = tempU.PDF(testCase.xvalues);
%                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-5,'RelTol',1e-5,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utConvolutionDisc


