classdef utChiSqNoncentral < utContinuous
    
    properties (ClassSetupParameter)
        parmdf     = struct( 'p1',1     , 'p2',2   , 'p4',4   , 'p8',8     , 'p16',16    , 'p32',32   );
        parmnoncen = struct( 'p_05',.05 , 'p_2',.2 , 'p_4',.4 , 'p_08',.08 , 'p1_6', 1.6 , 'p_25',0.25);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utChiSqNoncentral(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmdf,parmnoncen)
            % Computations specific to the ChiSqNoncentral distribution.
            testCase.Dist = ChiSqNoncentral(parmdf,parmnoncen);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % Adjust tolerances as appropriate for this distribution & parameters:

            testCase.EstParmCodes = 'fr';

            % testCase.Dist.UseSplineCDFOn(200);  % Approximation is quite poor due to fast rise from zero (case 1 at least)
            % testCase.Dist.UseSplinePDFOn(200);  % Approximation is quite poor due to high spike & long tail (case 1 at least)

            SetupXs(testCase,40,200);
            
            testCase.Expected.PDF = ncx2pdf( testCase.xvalues,parmdf,parmnoncen );
            testCase.Expected.CDF = ncx2cdf( testCase.xvalues,parmdf,parmnoncen );

            if parmdf > 1
                SetTolerances(testCase,0.0002);
                testCase.ParmEstAbsTol = 0.05;
                testCase.ParmEstRelTol = 0.01;
            else
                SetTolerances(testCase,0.002);
                testCase.CenMomentRelTol(3) = .005;
                testCase.RawMomentRelTol(4) = 0.01;
                testCase.KurtRelTol = 0.05;
                testCase.RawIntAbsTol(3) = 0.01;
                testCase.RawIntRelTol(4) = 0.015;
                testCase.RawIntRelTol(5) = 0.04;
                testCase.ParmEstAbsTol = 0.05;
                testCase.ParmEstRelTol = 0.05;
            end

%             testCase.Dist.SearchOptions.MaxFunEvals = 30000;
%             testCase.Dist.SearchOptions.MaxIter = 20000;
            testCase.Dist.SearchOptions.TolFun = 1e-8;
            testCase.Dist.SearchOptions.TolX = 1e-8;
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters.
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utChiSqNoncentral


