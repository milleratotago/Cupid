classdef utNormal < utContinuous
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'n100',-100 , 'n1',-1  , 'p0',0 , 'p5',5   , 'p25',25   , 'p250',250 );
        parmsigma = struct( 'p_01',.01  , 'p_1',.1 , 'p1',1 , 'p10',10 , 'p100',100 , 'p11',11   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utNormal(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmsigma)
            % Computations specific to the normal distribution.
            testCase.Dist = Normal(parmmu,parmsigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = (-3:.05:3)*parmsigma+parmmu;
            
            % Set up some X values for which MLE should return the true parameters:
            testCase.xMLE = -1:.1:1;
            StartingSD = std(testCase.xMLE,1);  % "1" indicates divisor in SD computation is N rather than N-1
            testCase.xMLE = testCase.xMLE * parmsigma / StartingSD + parmmu;
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.PDF = normpdf( (testCase.xvalues-parmmu)/parmsigma ) / parmsigma;
            testCase.Expected.CDF = normcdf( (testCase.xvalues-parmmu)/parmsigma );
            testCase.Expected.Mean = parmmu;
            testCase.Expected.Variance = parmsigma^2;
            testCase.Expected.SD = parmsigma;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            testCase.Expected.Kurtosis = 3;
            
            testCase.Expected.Median = parmmu;
            % testCase.Expected.Minimum =
            % testCase.Expected.Maximum =
            testCase.Expected.PctileSkew75 = 0;
            testCase.Expected.PctileSkew90 = 0;
            testCase.Expected.PctileSkew99 = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.RawSkewAbsTol = max(0.001, 0.01 * (abs(parmmu) + parmsigma));  % RawSkewness integral is problematic for large mu and/or sigma
            testCase.RawMomentAbsTol(4) = max(0.001, 0.01 * (abs(parmmu) + parmsigma));  % RawSkewness integral is problematic for large mu and/or sigma
            if parmsigma>1
                ExpansionFactor = parmsigma^0.25;   % Moments are less accurate when sigma is quite large
                testCase.CenMomentAbsTol = testCase.CenMomentAbsTol * ExpansionFactor;
                testCase.CenMomentRelTol = testCase.CenMomentRelTol * ExpansionFactor;
                testCase.ParmEstAbsTol = testCase.ParmEstAbsTol * ExpansionFactor;
                testCase.ParmEstRelTol = testCase.ParmEstRelTol * ExpansionFactor;
            end
            testCase.MLParmTolSE = nan;   % ML parameter estimation should be perfect

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utNormal


