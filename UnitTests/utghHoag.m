classdef utghHoag < utContinuous
    
    properties (ClassSetupParameter)
        parmA = struct( 'n10',-10    , 'n1',-1     , 'p0',0     , 'p5',5      , 'p25',25   , 'p2_5',2.5  );
        parmB = struct( 'p_5',.5     , 'p_1',.1    , 'p1',1     , 'p2',2      , 'p10',10   , 'p8'  ,8    );
        parmg = struct( 'p_1' ,.1    , 'n_02',-.02 , 'p_07',.07 , 'n_03',-.03 , 'p_08',.08 , 'n_04',-.04 );
        parmh = struct( 'p_001',.001 , 'p_1',.1    , 'p_2',.2  , 'p_05',.05  , 'p_07',.07 , 'p_08',.08  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utghHoag(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmA,parmB,parmg,parmh)
            % Computations specific to the ghHoag distribution.
            testCase.Dist = ghHoag(parmA,parmB,parmg,parmh);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            testCase.Dist.SplineXnPoints = 250;
            testCase.Dist.UseSplineTransXOn;
 
            SetupXs(testCase,40,2000);
            
            testCase.Expected.Median = parmA;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
%            testCase.RawMomentAbsTol(4) = max( [0.005, 0.001*abs(testCase.Dist.Mean), 0.001*abs(testCase.Dist.SD) ] );
             testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.001;
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.025;
             testCase.MLParmTolSE = 0.05;   % ML parameter estimation is not great
%            testCase.SkipEstAll = true;
             testCase.SkipEstMoment = true;
             testCase.SkipEstProbitAll = true;
             testCase.SkipEstPercentile = true;
             testCase.SkipEstChiSq = true;
             testCase.SkipEstPctBounds = true;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utghHoag


