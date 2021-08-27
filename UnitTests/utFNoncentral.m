classdef utFNoncentral < utContinuous
    
    properties (ClassSetupParameter)
        parm1DfNum   = struct( 'p1a',1    , 'p1b',1  , 'p1c',1    , 'p2',2      , 'p4',4      ,  'p8', 8     );
        parm2DfDenom = struct( 'p6',6     , 'p12',12 , 'p240',240 , 'p48', 48   , 'p96', 96   ,  'p320', 320 );
        parmnoncen   = struct( 'p_05',.05 , 'p_3',.3 , 'p1',1     , 'p_1',.1    , 'p_08',.08  ,  'p_01',.01  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utFNoncentral(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1DfNum,parm2DfDenom,parmnoncen)
            % Computations specific to the FNoncentral distribution.
            testCase.Dist = FNoncentral(parm1DfNum,parm2DfDenom,parmnoncen);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,1000);
            
            testCase.Expected.PDF = ncfpdf( testCase.xvalues,parm1DfNum,parm2DfDenom,parmnoncen );
            testCase.Expected.CDF = ncfcdf( testCase.xvalues,parm1DfNum,parm2DfDenom,parmnoncen );

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Not numerically very accurate
            testCase.KurtRelTol = .02;

            testCase.EstParmCodes = 'ffr';  % Not adjusting integer paramters
            testCase.ParmEstAbsTol(3) = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utFNoncentral

