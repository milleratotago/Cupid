classdef utF < utContinuous;
    
    properties (ClassSetupParameter)
        parm1DfNum   = struct( 'p1a',1 , 'p1b',1  , 'p1c',1    , 'p2',2      , 'p4',4       ,  'p8', 8   );
        parm2DfDenom = struct( 'p6',6  , 'p12',12 , 'p240',240 , 'p48', 48   ,  'p96', 96   ,  'p320', 320  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utF(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1DfNum,parm2DfDenom)
            % Computations specific to the F distribution.
            testCase.Dist = F(parm1DfNum,parm2DfDenom);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            % testCase.Dist.SearchOptions.TolFun = 1e-14;
            % testCase.Dist.SearchOptions.TolX = 1e-14;

            SetupXs(testCase,40,200);
            
            testCase.SkipMomentEst = true;  % Moments do not give much info about df's
            testCase.SkipMLEst = true;  % I can't find X values for which the MLE is the current dfs.
            % Set up some X values for which MLE should return (very close to) the true parameters:
            % Skipped because I can'F find X values for which the MLE is the current df.
%             npoints = 2000; % 20000;
%             testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.PDF = fpdf( testCase.xvalues, parm1DfNum,parm2DfDenom );
            testCase.Expected.CDF = fcdf( testCase.xvalues, parm1DfNum,parm2DfDenom );
            
            testCase.SkipMGFs = true;   % MGF does not exist for F distribution.

            testCase.HighestMoment = min(testCase.HighestMoment,floor((parm2DfDenom-0.5)/2));

            % Adjust tolerances as appropriate for this distribution & parameters:
            if parm1DfNum > 1
                SetTolerances(testCase,0.025);
            else
                SetTolerances(testCase,0.05);
            end
            testCase.CenMomentRelTol(3) = .05;   % Numerical integration is not very accurate.
            testCase.CenMomentRelTol(4) = .05;
            testCase.RawMomentRelTol(3) = .05;
            testCase.RawMomentRelTol(4) = .05;
            testCase.KurtRelTol = 0.1;
            if testCase.HighestMoment+1>=3
                testCase.RawIntAbsTol(3) = 0.01;
            end
            if testCase.HighestMoment+1>=4
                testCase.RawIntRelTol(4) = 0.015;
            end
            if testCase.HighestMoment+1>=5
                testCase.RawIntRelTol(5) = 0.04;
            end
%             testCase.ParmEstAbsTol = 1;  % OK if df estimates are off by 1.
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters.
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utt


