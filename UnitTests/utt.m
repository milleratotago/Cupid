classdef utt < utContinuous;
    
    properties (ClassSetupParameter)
        parmDF = struct( 'p6',6 , 'p12',12 , 'p24',24 ,  'p48', 48   ,  'p96', 96   ,  'p320', 320  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utt(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmDF)
            % Computations specific to the t distribution.
            testCase.Dist = t(parmDF);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            testCase.Dist.SearchOptions.MaxIter = 20000;
%             testCase.Dist.SearchOptions.TolFun = 1e-14;
%             testCase.Dist.SearchOptions.TolX = 1e-14;

            SetupXs(testCase,40,200);
            
            % Set up some X values for which MLE should return (very close to) the true parameters:
            % Skipped because I can't find X values for which the MLE is the current df.
            % npoints = 2000; % 20000;
            % testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.PDF = tpdf( testCase.xvalues, parmDF );
            testCase.Expected.CDF = tcdf( testCase.xvalues, parmDF );
            testCase.Expected.Mean = 0;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            
            testCase.MGFh = 1.0E-4; % Too large or too small produces numerical errors.

            % Adjust tolerances as appropriate for this distribution & parameters:
            if parmDF > 1
                SetTolerances(testCase,0.002);
            else
                SetTolerances(testCase,0.01);
                testCase.CenMomentRelTol(3) = .005;
                testCase.RawMomentRelTol(4) = 0.01;
                testCase.KurtRelTol = 0.05;
                testCase.RawIntAbsTol(3) = 0.01;
                testCase.RawIntRelTol(4) = 0.015;
                testCase.RawIntRelTol(5) = 0.04;
            end
            testCase.SkipMLEst = true;  % I can't find X values for which the MLE is the current df.
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters.
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utt


