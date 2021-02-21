classdef utGeneric < matlab.unittest.TestCase
    
    % An abstract class for distribution checking.
    %    Checks specific known values computed externally if they are known.
    %    Performs a variety of internal consistency tests that should be passed by _every_ probability distribution.
    
    properties (Abstract)
        % Use something Abstract so that matlab.unittest.TestSuite.fromFolder
        % will NOT attempt to define tests from this file.
        Dummy1
    end
    
    properties
        Dist           % The probability distribution being tested.
        xvalues        % Value at which the distribution is to be evaluated/checked; this must be set by each child distribution.
        Computed       % A structure to hold various values computed from the distribution.
        Expected       % A structure to hold various values that are known "a priori" (i.e., from sources outside Cupid).
        
        % **************** Parameters used to control testing:
        
        HighestMoment  % Moments are computed from 0 to the highest one.
        
        MGFh            % A small constant used in computing derivatives of MGF
        
        % Parameters controlling the chi-square test for random numbers fitting the distribution:
        ChiSqNTries    % Number of tries to get random numbers passing chi-square test
        ChiSqCriticalp % p value of the chi-square test must be at least this large to pass.
        ChiSqNRands    % Number of random numbers to generate
        ChiSqNBins     % N of bins to use in chi-square test of random number generator
        
        % Variables used during parameter estimation:
        EstParmCodes   % Parmcodes indicating whether each parameter is real or integer or fixed
        
        xMLE           % X values to be used in MLE parameter estimation.  Must be set for each distribution.
        ChiSqEstNBins  % N of bins used in ChiSqEst testing.
        
        % Parameters that can be set to control probit estimation:
        ProbitNBins    % N of bins to be considered.
        ProbitNPerBin  % N of trials per bin to be considered.
        ProbitmAFC     % N of alternatives in mAFC task to be modeled.
        % These parameters are also used in probit estimation, but they are calculated automatically.
        ProbitBinMax, ProbitBinCDFs, ProbitNTrials
        
        % Absolute & relative tolerances for cross-checking computed/expected values with verifyEqual (& CheckIfExpectedKnown):
        CDFAbsTol
        CDFRelTol
        CenIntAbsTol      % Vector with HighestMoment+1 elements
        CenIntRelTol      % Vector with HighestMoment+1 elements
        CenMomentAbsTol   % Vector with HighestMoment+1 elements
        CenMomentRelTol   % Vector with HighestMoment+1 elements
        HazAbsTol
        HazRelTol
        KurtAbsTol
        KurtRelTol
        MeanAbsTol
        MeanRelTol
        MGFMom2AbsTol
        MGFMom2RelTol
        MGFXAbsTol
        MGFXRelTol
        MLParmTolSE   % If this number is specified, ML estimates are accepted if they are within
        % SE*MLParmTolSE of the true values (or if they are within ParmEstAbsTol or Rel).
        ParmEstAbsTol     % Vectors with NDistParms elements.  Parameter estimation processes
        ParmEstRelTol     % should return parameters within this close to their original values.
        PctileSkewAbsTol
        PctileSkewRelTol
        PDFAbsTol
        PDFRelTol
        RawIntAbsTol      % Vector with HighestMoment+1 elements
        RawIntRelTol      % Vector with HighestMoment+1 elements
        RawMomentAbsTol   % Vector with HighestMoment+1 elements
        RawMomentRelTol   % Vector with HighestMoment+1 elements
        RawSkewAbsTol
        RawSkewRelTol
        RelSkewAbsTol
        RelSkewRelTol
        VarAbsTol
        VarRelTol
        XAbsTol
        XRelTol
        MLErrTol      % Estimation may be allowed to increase error slightly.
        PDF0AbsTol, CDF01AbsTol   % For checking that PDF is zero outside bounds and CDF is 0/1
        
        % Minimum and maximum allowable random numbers:
        MinRand, MaxRand
        
        % Flags to skip certain tests
        SkipMGFs
        SkipEstML
        SkipEstPercentile
        SkipEstChiSq
        SkipEstPctBounds
        SkipEstMoment
        SkipEstAll
        SkipEstProbitAll
    end
    
    methods
        
        function testCase=utGeneric(varargin)  % Constructor
            testCase=testCase@matlab.unittest.TestCase;

            % Default property initializations:
            testCase.HighestMoment = 4;
            
            % Chapter 5.7 of Numerical Recipes in C (http:%www.nrbook.com/a/bookcpdf/c5-7.pdf)
            % recommends a different (smaller) choice for MGFh, but I find their recommendation
            % is a bit too small for many distributions.
            testCase.MGFh = 1.0E-6; % Too large or too small produces numerical errors. This value works quite well for normal.
            
            testCase.ChiSqNRands = 10000;
            testCase.ChiSqNBins = 100;
            testCase.ChiSqNTries = 3;
            testCase.ChiSqCriticalp = .05;
            
            testCase.ChiSqEstNBins = 50;
            
            testCase.ProbitNBins = 10;
            testCase.ProbitNPerBin = 50;
            testCase.ProbitmAFC = 2;
            
            testCase.SkipMGFs = false;

            testCase.SkipEstML = false;
            testCase.SkipEstPercentile = false;
            testCase.SkipEstChiSq = false;
            testCase.SkipEstPctBounds = false;
            testCase.SkipEstMoment = false;
            testCase.SkipEstAll = false;
            testCase.SkipEstProbitAll = false;
        end
        
        function SetTolerances(testCase,DefaultTolerance)
            % Set default absolute & relative tolerances for checking values with verifyEqual (& CheckIfExpectedKnown):
            % DefaultTolerance = sqrt(eps);   % This would be a very, very, very tough tolerance!
            %             DefaultTolerance = 0.0001;      % Perhaps this will be appropriate in many cases.
            testCase.CDFAbsTol = DefaultTolerance;
            testCase.CDFRelTol = DefaultTolerance;
            testCase.CenIntAbsTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.CenIntRelTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.CenMomentAbsTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.CenMomentRelTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.HazAbsTol = DefaultTolerance;
            testCase.HazRelTol = DefaultTolerance;
            testCase.KurtAbsTol = DefaultTolerance;
            testCase.KurtRelTol = DefaultTolerance;
            testCase.MeanAbsTol = DefaultTolerance;
            testCase.MeanRelTol = DefaultTolerance;
            testCase.MGFMom2AbsTol = DefaultTolerance;
            testCase.MGFMom2RelTol = DefaultTolerance * 3;  % RELAX
            testCase.MGFXAbsTol = DefaultTolerance;
            testCase.MGFXRelTol = DefaultTolerance;
            % Note that ParmEstTol cannot be defined until after Dist has been constructed.
            testCase.EstParmCodes = testCase.EstParmCodes;
            testCase.MLParmTolSE = 0.1;  % Check relative to SE of MLE's by default
            testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * DefaultTolerance;
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * DefaultTolerance;
            testCase.PctileSkewAbsTol = DefaultTolerance;
            testCase.PctileSkewRelTol = DefaultTolerance;
            testCase.PDFAbsTol = DefaultTolerance;
            testCase.PDFRelTol = DefaultTolerance;
            testCase.RawIntAbsTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.RawIntRelTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.RawMomentAbsTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.RawMomentRelTol = ones(1,testCase.HighestMoment+1) * DefaultTolerance;
            testCase.RawSkewAbsTol = DefaultTolerance;
            testCase.RawSkewRelTol = DefaultTolerance;
            testCase.RelSkewAbsTol = DefaultTolerance;
            testCase.RelSkewRelTol = DefaultTolerance;
            testCase.VarAbsTol = DefaultTolerance;
            testCase.VarRelTol = DefaultTolerance;
            testCase.XAbsTol = DefaultTolerance;
            testCase.XRelTol = DefaultTolerance;
            testCase.MLErrTol = 0;
            testCase.PDF0AbsTol = 0.0001;
            testCase.CDF01AbsTol = 0.0001;
            
            testCase.MinRand = testCase.Dist.LowerBound;
            testCase.MaxRand = testCase.Dist.UpperBound;

            % *** CONTROL SPEED HERE ***
            testCase.Dist.SearchOptions.TolX = 1e-3;  % Control accuracy in parameter searching for estimation (ML, probit, etc)
            
        end
        
        
        function CheckIfExpectedKnown(testCase,s,AbsTol,RelTol)
            if isfield(testCase.Expected,s)
                testCase.verifyEqual(testCase.Computed.(s),testCase.Expected.(s),'AbsTol',AbsTol,'RelTol',RelTol,['Computed vs Expected ' s]);
            end
        end
        
        function NCorrect = MakeProbitmAFCData(testCase)
            % The data to be fit are a set of bins defined by the original distribution,
            % and the value of CDF for each of those bins, adjusted for guessing.
            GuessProb = 1 / testCase.ProbitmAFC;
            NCorrect = testCase.ProbitNTrials .* (GuessProb + (1 - GuessProb) * testCase.ProbitBinCDFs);
        end
        
        function utGenericMethodSetup(testCase)
            
            if numel(testCase.EstParmCodes)==0
                testCase.EstParmCodes = testCase.Dist.DefaultParmCodes;
            end
            
            global WantPlots;
            if WantPlots
                testCase.Dist.PlotDens;
                drawnow;
            end
            
            testCase.Computed.PDF = testCase.Dist.PDF(testCase.xvalues);
            testCase.Computed.CDF = testCase.Dist.CDF(testCase.xvalues);
            testCase.Computed.Hazard = testCase.Dist.Hazard(testCase.xvalues);
            testCase.Computed.CDF(testCase.Computed.CDF<0) = 0;
            testCase.Computed.CDF(testCase.Computed.CDF>1) = 1;
            testCase.Computed.InverseCDF = testCase.Dist.InverseCDF(testCase.Computed.CDF);
%            testCase.Computed.InverseCDF = testCase.Dist.InverseCDF(targetcdf);
            testCase.Expected.InverseCDF = testCase.xvalues;
            
            % Compute moments:
            testCase.Computed.RawMoment = zeros(1,testCase.HighestMoment+1);
            testCase.Computed.CenMoment = zeros(1,testCase.HighestMoment+1);
            for iPower=0:testCase.HighestMoment
                testCase.Computed.RawMoment(iPower+1) = testCase.Dist.RawMoment(iPower);
                testCase.Computed.CenMoment(iPower+1) = testCase.Dist.CenMoment(iPower);
            end
            
            % Compute functions of moments:
            testCase.Computed.Mean = testCase.Dist.Mean;
            testCase.Computed.Variance = testCase.Dist.Variance;
            testCase.Computed.SD = testCase.Dist.SD;
            testCase.Computed.RawSkewness = testCase.Dist.RawSkewness;
            testCase.Computed.RelSkewness = testCase.Dist.RelSkewness;
            testCase.Computed.Kurtosis = testCase.Dist.Kurtosis;
            
            testCase.Computed.Median = testCase.Dist.Median;
            testCase.Computed.Minimum = testCase.Dist.Minimum;
            testCase.Computed.Maximum = testCase.Dist.Maximum;
            testCase.Computed.PctileSkew75 = testCase.Dist.PctileSkew(0.75);
            testCase.Computed.PctileSkew90 = testCase.Dist.PctileSkew(0.90);
            testCase.Computed.PctileSkew99 = testCase.Dist.PctileSkew(0.99);
            
            % Make bins for probit testing:
            testCase.ProbitBinMax = testCase.Dist.MakeBinSet(1/testCase.ProbitNBins);
            testCase.ProbitBinMax = testCase.ProbitBinMax(2:end-1);  % omit the 2 most extreme bins to avoid edge effects
            testCase.ProbitBinCDFs = testCase.Dist.CDF(testCase.ProbitBinMax);
            testCase.ProbitNTrials = ones(size(testCase.ProbitBinCDFs))*testCase.ProbitNPerBin;
            
        end
        
        function tf = TooFewForProbit(testCase,Complain)       
            tf = (testCase.Dist.DistType == 'd') && (testCase.ProbitNBins >= testCase.Dist.NValues);
            if tf && Complain
                 warning('Not enough discrete values to do probit estimation.');
            end
        end

    end  % methods
    
    methods (Test, ParameterCombination='sequential')
        
        function PDFcheck(testCase)
            % Make sure the PDF can be computed & is always positive.
            testCase.verifyGreaterThanOrEqual(testCase.Computed.PDF,0,'Found negative PDF values');
            CheckIfExpectedKnown(testCase,'PDF',testCase.PDFAbsTol,testCase.PDFRelTol);
            % Make sure the PDF is zero for values out of bounds.
            OutOfRange = [testCase.Dist.LowerBound-(1:3) testCase.Dist.UpperBound+(1:3)];
            testCase.verifyEqual(testCase.Dist.PDF(OutOfRange),zeros(size(OutOfRange)),'AbsTol',testCase.PDF0AbsTol,'Found non-zero PDF values for out-of-range X values.');
        end
        
        function CDFcheck(testCase)
            % Make sure the CDF can be computed & is always in the range 0-1.
            testCase.verifyGreaterThanOrEqual(testCase.Computed.CDF,0,'Found negative CDF values');
            testCase.verifyLessThanOrEqual(testCase.Computed.CDF,1,'Found CDF values greater than 1');
            CheckIfExpectedKnown(testCase,'CDF',testCase.CDFAbsTol,testCase.CDFRelTol);
            % Make sure the CDF is zero for values below the lower bound and is one for values above the upper bound.
            OutOfRange = [testCase.Dist.LowerBound-(1:3)];
            testCase.verifyEqual(testCase.Dist.CDF(OutOfRange),zeros(size(OutOfRange)),'AbsTol',testCase.CDF01AbsTol,'Found non-zero CDF values for X values below lower bound.');
            OutOfRange = [testCase.Dist.UpperBound+(1:3)];
            testCase.verifyEqual(testCase.Dist.CDF(OutOfRange),ones(size(OutOfRange)),'AbsTol',testCase.CDF01AbsTol,'Found non-one CDF values for X values above upper bound.');
        end
        
        function Hazardcheck(testCase)
            % Make sure the hazard function can be computed & is always positive.
            Useable = testCase.Computed.CDF < 1;
            if numel(testCase.Computed.Hazard(Useable))>0  % may be zero for constant distributions
                testCase.verifyGreaterThanOrEqual(testCase.Computed.Hazard(Useable),0,'Found negative Hazard function values');
            end
            % Check that the Hazard function values are consistent with the PDF/CDF values:
            HazardFromPDFCDF = testCase.Computed.PDF(Useable) ./ (1 - testCase.Computed.CDF(Useable));
            testCase.verifyEqual(testCase.Computed.Hazard(Useable),HazardFromPDFCDF,'AbsTol',testCase.HazAbsTol,'RelTol',testCase.HazRelTol,'Hazard function values are not consistent with PDF/CDF values within HazAbsTol/HazRelTol.');
            CheckIfExpectedKnown(testCase,'Hazard',testCase.HazAbsTol,testCase.HazRelTol);
        end
        
        function InverseCDFcheck(testCase)
            % Make sure the CDF can be computed & is always in the range 0-1.
            testCase.verifyEqual(testCase.Computed.InverseCDF,testCase.xvalues,'AbsTol',testCase.XAbsTol,'RelTol',testCase.XRelTol,'InverseCDF values do not match original xvalues within XAbsTol/XRelTol.');
            % CheckIfExpectedKnown(testCase,'InverseCDF',testCase.XAbsTol,testCase.XRelTol);  % This test would be redundant because Expected was simply set to xvalues
        end
        
        function Moments(testCase)
            
            % Check for consistency of moments & functions of moments:
            testCase.verifyEqual(testCase.Computed.RawMoment(1),1,'AbsTol',testCase.RawMomentAbsTol(1),'RelTol',testCase.RawMomentRelTol(1),'RawMoment integral ^0 unequal to 1 within RawMomentAbsTol(1).');
            testCase.verifyEqual(testCase.Computed.CenMoment(1),1,'AbsTol',testCase.CenMomentAbsTol(1),'RelTol',testCase.CenMomentRelTol(1),'CenMoment integral ^0 unequal to 1 within CenMomentAbsTol(1).');
            if testCase.HighestMoment>=1
                testCase.verifyEqual(testCase.Computed.RawMoment(2),testCase.Computed.Mean,'AbsTol',testCase.RawMomentAbsTol(2),'RelTol',testCase.RawMomentRelTol(2),'RawMoment integral ^1 unequal to mean within RawMomentAbsTol(2).');
                testCase.verifyEqual(testCase.Computed.CenMoment(2),0,'AbsTol',testCase.CenMomentAbsTol(2),'RelTol',testCase.CenMomentRelTol(2),'CenMoment integral ^1 unequal to 0 within CenMomentAbsTol(2).');
            end
            if testCase.HighestMoment>=2
                VarFromRaw = testCase.Computed.RawMoment(3) - testCase.Computed.RawMoment(2)^2;
                testCase.verifyEqual(VarFromRaw,testCase.Computed.Variance,'AbsTol',testCase.RawMomentAbsTol(3),'RelTol',testCase.RawMomentRelTol(3),'RawMoment integrals ^1 & ^2 do not give computed Variance within RawMomentAbsTol(3).');
                testCase.verifyEqual(testCase.Computed.CenMoment(3),testCase.Computed.Variance,'AbsTol',testCase.CenMomentAbsTol(3),'RelTol',testCase.CenMomentRelTol(3),'CenMoment integral ^2 unequal to Variance within CenMomentAbsTol(3).');
            end
            if testCase.HighestMoment>=3
                SkewFromMom = testCase.Computed.RawMoment(4) - 3*testCase.Computed.RawMoment(2)*testCase.Computed.CenMoment(3) - testCase.Computed.RawMoment(2)^3;
                if SkewFromMom >= 0
                    SkewFromMom = SkewFromMom^(1/3);
                else
                    SkewFromMom = -(-SkewFromMom)^(1/3);
                end
                testCase.verifyEqual(SkewFromMom,testCase.Computed.RawSkewness,'AbsTol',testCase.RawMomentAbsTol(4),'RelTol',testCase.RawMomentRelTol(4),'RawSkewness is not consistent with RawMoments ^1 to ^3 within RawMomentAbsTol(4).');
                testCase.verifyEqual(testCase.Computed.RawSkewness^3,testCase.Computed.CenMoment(4),'AbsTol',testCase.CenMomentAbsTol(4),'RelTol',testCase.CenMomentRelTol(4),'RawSkewness^3 does not equal CenMoment^3 within CenMomentAbsTol(4).');
            end
            if testCase.HighestMoment>=4
                testCase.verifyEqual(testCase.Computed.Kurtosis,testCase.Computed.CenMoment(5)/testCase.Computed.CenMoment(3)^2,'AbsTol',testCase.KurtAbsTol,'RelTol',testCase.KurtRelTol,'Kurtosis does not equal CenMoment^4/CenMoment^2 within KurtAbsTol/KurtRelTol.');
            end
            
            % Check that computed values match any values that are expected a priori.
            CheckIfExpectedKnown(testCase,'Mean',testCase.MeanAbsTol,testCase.MeanRelTol);
            CheckIfExpectedKnown(testCase,'Variance',testCase.VarAbsTol,testCase.VarRelTol);
            CheckIfExpectedKnown(testCase,'SD',sqrt(testCase.VarAbsTol),sqrt(testCase.VarRelTol));
            CheckIfExpectedKnown(testCase,'RawSkewness',testCase.RawSkewAbsTol,testCase.RawSkewRelTol);
            CheckIfExpectedKnown(testCase,'RelSkewness',testCase.RelSkewAbsTol,testCase.RelSkewRelTol);
            CheckIfExpectedKnown(testCase,'Kurtosis',testCase.KurtAbsTol,testCase.KurtRelTol);
            CheckIfExpectedKnown(testCase,'RawMoment',testCase.RawMomentAbsTol,testCase.RawMomentRelTol);
            CheckIfExpectedKnown(testCase,'CenMoment',testCase.CenMomentAbsTol,testCase.CenMomentRelTol);
            
            CheckIfExpectedKnown(testCase,'Median',testCase.XAbsTol,testCase.XRelTol);
            CheckIfExpectedKnown(testCase,'Minimum',testCase.XAbsTol,testCase.XRelTol);
            CheckIfExpectedKnown(testCase,'Maximum',testCase.XAbsTol,testCase.XRelTol);
            CheckIfExpectedKnown(testCase,'PctileSkew75',testCase.PctileSkewAbsTol,testCase.PctileSkewRelTol);
            CheckIfExpectedKnown(testCase,'PctileSkew90',testCase.PctileSkewAbsTol,testCase.PctileSkewRelTol);
            CheckIfExpectedKnown(testCase,'PctileSkew99',testCase.PctileSkewAbsTol,testCase.PctileSkewRelTol);
            
        end
        
        function MGFs(testCase)
            if testCase.SkipMGFs
                return
            end
            MGF0 = testCase.Dist.MGF(0);  % Should be 1 but there can be numerical problems.
            if isnan(MGF0)
                MGF0 = 1;
            end
            MGFPlusH = testCase.Dist.MGF(testCase.MGFh);
            MGFMinusH = testCase.Dist.MGF(-testCase.MGFh);
            MeanByMGF = (MGFPlusH - MGFMinusH) / (2*testCase.MGFh);
            Mom2ByMGF = (MGFPlusH + MGFMinusH - 2*MGF0) / testCase.MGFh^2;
            if isfield(testCase.Expected,'Mean')
                ExpectedMean = testCase.Expected.Mean;
            else
                ExpectedMean = testCase.Computed.Mean;
            end
            testCase.verifyEqual(ExpectedMean,MeanByMGF,'AbsTol',testCase.MGFXAbsTol,'RelTol',testCase.MGFXRelTol,'Mean does not match that estimated from MGF within MGFXAbsTol/MGFXRelTol.');
            if isfield(testCase.Expected,'Variance')
                ExpectedVariance = testCase.Expected.Variance;
            else
                ExpectedVariance = testCase.Computed.Variance;
            end
            testCase.verifyEqual(ExpectedVariance+ExpectedMean^2,Mom2ByMGF,'AbsTol',testCase.MGFMom2AbsTol,'RelTol',testCase.MGFMom2RelTol,'Moment2 does not match that estimated from MGF within MGFMom2AbsTol/MGFMom2RelTol.');
        end  % MGFs
        
        function CheckRandom(testCase)
            % Check that samples of random numbers are consistent with the PDF of the
            % distribution by a chi-square test.
            % Skip if it is a discrete distribution with only 1 value.
            if (testCase.Dist.DistType=='d') && (testCase.Dist.NValues==1)
                return;
            end
            PassedChisqTest = false;
            NTries = 0;
            if testCase.Dist.DistType=='c'
                [BinMax, BinProb] = testCase.Dist.MakeBinSet(1/testCase.ChiSqNBins);
            end
% NWJEFF: Including this caused problems with LinearTrans(Binomial(53,0.44),-2.3,55)--obschisq histcounts found too many observations in smallish bins.
%         Omitting it caused problems with all utBinomial, I guess because histcounts found not enough observations in smallish bins
% I think I fixed this by adding ObsChiSq to dDiscrete
%             if testCase.Dist.DistType == 'd'  
%                 BinMax = BinMax + eps(BinMax);
%             end
            while ~PassedChisqTest && (NTries < testCase.ChiSqNTries)
                NTries = NTries + 1;
                RandVals = testCase.Dist.Random(testCase.ChiSqNRands,1);
                testCase.verifyGreaterThanOrEqual(RandVals,testCase.MinRand,'Generated random number(s) below the lower bound.');
                testCase.verifyLessThanOrEqual(RandVals,testCase.MaxRand,'Generated random number(s) greater than the upper bound.');
                if testCase.Dist.DistType=='c'
                    [~, obschisqp] = obschisq(RandVals,BinMax,BinProb);
                elseif testCase.Dist.DistType=='d'
                    [~, obschisqp] = testCase.Dist.ObsChiSq(RandVals);
                else
                    error('CheckRandom has not been implemented for mixed discrete/continuous distributions.');
                end
                PassedChisqTest = obschisqp > testCase.ChiSqCriticalp;
            end
            testCase.verifyTrue(PassedChisqTest,'Random numbers failed the chi-square test.');
        end % CheckRandom
        
        function CheckIntegrals(testCase)
            RawInt = zeros(1,testCase.HighestMoment+1);
            CenInt = zeros(1,testCase.HighestMoment+1);
            for iPwr = 0:testCase.HighestMoment
                Ip1 = iPwr+1;
                RawInt(Ip1) = testCase.Dist.IntegralXToNxPDF(testCase.Dist.LowerBound,testCase.Dist.UpperBound,iPwr);
                CenInt(Ip1) = testCase.Dist.IntegralX_CToNxPDF(testCase.Dist.LowerBound,testCase.Dist.UpperBound,testCase.Computed.Mean,iPwr);
            end
            testCase.verifyEqual(RawInt,testCase.Computed.RawMoment,'AbsTol',testCase.RawIntAbsTol,'RelTol',testCase.RawIntRelTol,'Raw integrals do not match moments within RawIntAbsTol/RawIntRelTol.');
            testCase.verifyEqual(CenInt,testCase.Computed.CenMoment,'AbsTol',testCase.CenIntAbsTol,'RelTol',testCase.CenIntRelTol,'Central integrals do not match moments within CenIntAbsTol/CenIntRelTol.');
        end  % CheckIntegrals
        
        function MLEstTest(testCase)
            % The data to be fit are the original testCase.xMLE values.

            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstML
                return
            end
            
            HoldWarn = testCase.Dist.SkipImpossibleWarn;
            testCase.Dist.SkipImpossibleWarn = true;  % Turn off warnings about impossible data values
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            %             t = size(testCase.xMLE)
            %             s=testCase.Dist.StringName
            ErrorWithTrueParms = -testCase.Dist.LnLikelihood(testCase.xMLE);
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = -testCase.Dist.LnLikelihood(testCase.xMLE);
            testCase.Dist.EstML(testCase.xMLE,testCase.EstParmCodes);
            EndError = -testCase.Dist.LnLikelihood(testCase.xMLE);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'MLE estimation increased error.');
            
            % Check whether the original parameter values were approximately recovered as they should have been if xvalues were chosen correctly.
            if isnan(testCase.MLParmTolSE)
                MaxAbsTol = testCase.ParmEstAbsTol;
            else
                FracSEs = testCase.Dist.MLSE(testCase.xMLE,testCase.EstParmCodes)*testCase.MLParmTolSE;
                MaxAbsTol = ones(size(testCase.EstParmCodes)).*testCase.ParmEstAbsTol;
                % FracSEs may be a different size if some solutions have imaginary components.
                if isreal(FracSEs) && sum(isinf(FracSEs))==0
                    MaxAbsTol = max(MaxAbsTol,FracSEs);
                end
            end
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',MaxAbsTol,'RelTol',testCase.ParmEstRelTol,'MLE estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            % The next line cannot be used because I do not know how (in general) to choose x's such that the true parameters
            % necessarily have the minimum error.
            %           testCase.verifyGreaterThanOrEqual(EndError,ErrorWithTrueParms-testCase.MLErrTol,'xMLE values do not produce minimal error at the true parameter values.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
            testCase.Dist.SkipImpossibleWarn = HoldWarn;
        end
        
        function MomEstTest(testCase)
            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstMoment
                return
            end
            
            Adjustable = (testCase.EstParmCodes == 'r') | (testCase.EstParmCodes == 'i');
            NFreeParms = sum(Adjustable);
            
            % The data to be fit are the true moments.
            % Always take at least 2 moments because t distribution has 1 parameter yet mean always 0.  Doesn't create problems for exponential.
            NMom = min(numel(testCase.Computed.CenMoment)-1,max(2,NFreeParms)); % NWJEFF: max (1+obj.NConstrainedMoments, NFree) ???
            ObsMoments = testCase.Computed.CenMoment(2:NMom+1);
            ObsMoments(1) = testCase.Computed.RawMoment(2);
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = testCase.Dist.MomentError(ObsMoments);
            testCase.Dist.EstMom(ObsMoments,testCase.EstParmCodes);
            EndError = testCase.Dist.MomentError(ObsMoments);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'Moment estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'Moment estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms); % Reinstate distribution parameters at the end of each fit.
        end
        
        function PercentileEstTest(testCase)
            % The data to be fit are the original testCase.xvalues & their CDF values.
            
            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstPercentile
                return
            end
            
            ObsPercentiles = testCase.xvalues;
            TargetCDF = testCase.Computed.CDF;
            %           fprintf('%8.3f %8.3f\n',[ObsPercentiles',TargetCDF']');  % testing
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = testCase.Dist.PercentileError(ObsPercentiles,TargetCDF);
            testCase.Dist.EstPctile(ObsPercentiles,TargetCDF,testCase.EstParmCodes);
            EndError = testCase.Dist.PercentileError(ObsPercentiles,TargetCDF);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'Percentile estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            % With discrete distributions, parameter values may not be recovered very well because of graininess in the X value
            % and jumpiness in the CDF function.  To compensate for that, we accept the estimation if the new CDF values are
            % very close to correct, even if the parameters are different.
            if (testCase.Dist.DistType == 'd') && (EndError < 0.01)
                CheckParms = HoldParms;
            else
                CheckParms = testCase.Dist.ParmValues;
            end
            testCase.verifyEqual(CheckParms,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'Percentile estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');

            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end
        
        function PctBoundsEstTest(testCase)
            % The data to be fit are two percentile points of the true distribution.
            
            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstPctBounds
                return
            end
            
            LowerP = 0.2;
            UpperP = 0.8;
            LowerX = testCase.Dist.InverseCDF(LowerP);
            UpperX = testCase.Dist.InverseCDF(UpperP);
            TargetpDiff = UpperP - LowerP;
            %           fprintf('%8.3f %8.3f\n',[ObsPercentiles',TargetCDF']');  % testing
            if LowerX >= UpperX
                warning('Cannot perform PctBoundsEstTest for this distribution because LowerX >= UpperX');
                return;
            end
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = abs(TargetpDiff - testCase.Dist.CDF(UpperX) + testCase.Dist.CDF(LowerX));
            testCase.Dist.EstPctBounds(LowerX,UpperX,TargetpDiff,testCase.EstParmCodes);
            EndError = abs(TargetpDiff - testCase.Dist.CDF(UpperX) + testCase.Dist.CDF(LowerX));
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'Pct bounds estimation increased error.');
            
            % There is no point in checking whether the original parameter values were recovered, because they need not have been.
            % With most distributions, various different parameter combinations yield the target percentage between LowerX and UpperX,
            % so there is no reason to expect that the original parameter values will be recovered.
            % testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'Pct bounds estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');

            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end
        
        function ChiSqEstTest(testCase)
            % The data to be fit are a set of bins determined by the original distribution, and the associated bin probabilities.
            
            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstChiSq
                return
            end
            
            BinMax = testCase.Dist.MakeBinSet(1/testCase.ChiSqEstNBins);
            BinCDFs = testCase.Dist.CDF(BinMax);
            BinProbs = diff([0 BinCDFs]);
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = testCase.Dist.GofFChiSq(BinMax,BinProbs);
            testCase.Dist.EstChiSq(BinMax,BinProbs,testCase.EstParmCodes);
            EndError = testCase.Dist.GofFChiSq(BinMax,BinProbs);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'ChiSq estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'ChiSq estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end

        function ProbitYNMaxLikEstTest(testCase)
            % The data to be fit are a set of bins defined by the original distribution, and the value of 1-CDF for each of those bins.
            
            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstProbitAll || TooFewForProbit(testCase,true)
                return
            end

            NGreater = testCase.ProbitNTrials .* testCase.ProbitBinCDFs;
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = -testCase.Dist.YNProbitLnLikelihood(testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.Dist.EstProbitYNML(testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater,testCase.EstParmCodes);
            EndError = -testCase.Dist.YNProbitLnLikelihood(testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'ProbitYNMaxLikEst estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'ProbitYNMaxLikEst estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end % ProbitYNMaxLikEstTest
        
        function ProbitYNChiSqEstTest(testCase)
            % The data to be fit are a set of bins defined by the original distribution, and the value of 1-CDF for each of those bins.
            
            global GlobalSkipEstAll;
            if GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstProbitAll || TooFewForProbit(testCase,false)
                return
            end
            
            NGreater = testCase.ProbitNTrials .* testCase.ProbitBinCDFs;
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = testCase.Dist.YNProbitChiSq(testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.Dist.EstProbitYNChiSq(testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater,testCase.EstParmCodes);
            EndError = testCase.Dist.YNProbitChiSq(testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'ProbitYNChiSqEst estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'ProbitYNChiSqEst estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end % ProbitYNChiSqEstTest
        
        function ProbitmAFCMaxLikEstTest(testCase)
            % The data to be fit are a set of bins defined by the original distribution, and the value of 1-CDF for each of those bins.
            
            global GlobalSkipEstAll;
            if  GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstProbitAll || TooFewForProbit(testCase,false)
                return
            end
            
            NGreater = MakeProbitmAFCData(testCase);
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = -testCase.Dist.mAFCProbitLnLikelihood(testCase.ProbitmAFC,testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.Dist.EstProbitmAFCML(testCase.ProbitmAFC,testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater,testCase.EstParmCodes);
            EndError = -testCase.Dist.mAFCProbitLnLikelihood(testCase.ProbitmAFC,testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'ProbitmAFCMaxLikEst estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'ProbitmAFCMaxLikEst estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end % ProbitmAFCMaxLikEstTest
        
        function ProbitmAFCChiSqEstTest(testCase)
            % The data to be fit are a set of bins defined by the original distribution, and the value of 1-CDF for each of those bins.
            
            global GlobalSkipEstAll;
            if  GlobalSkipEstAll || testCase.SkipEstAll || testCase.SkipEstProbitAll || TooFewForProbit(testCase,false)
                return
            end
            
            NGreater = MakeProbitmAFCData(testCase);
            
            % Save current parameters & then perturb them slightly, to see whether estimation process recovers the original parms:
            HoldParms = testCase.Dist.ParmValues;
            testCase.Dist.PerturbParms(testCase.EstParmCodes);
            
            % Check error with the perturbed parameters, try to find the new ones, and make sure error was reduced.
            StartError = testCase.Dist.mAFCProbitChiSq(testCase.ProbitmAFC,testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.Dist.EstProbitmAFCChiSq(testCase.ProbitmAFC,testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater,testCase.EstParmCodes);
            EndError = testCase.Dist.mAFCProbitChiSq(testCase.ProbitmAFC,testCase.ProbitBinMax,testCase.ProbitNTrials,NGreater);
            testCase.verifyGreaterThanOrEqual(StartError,EndError-testCase.MLErrTol,'ProbitmAFCChiSqEst estimation increased error.');
            
            % Check whether the original parameter values were recovered as they should have been.
            testCase.verifyEqual(testCase.Dist.ParmValues,HoldParms,'AbsTol',testCase.ParmEstAbsTol,'RelTol',testCase.ParmEstRelTol,'ProbitmAFCChiSqEst estimation did not recover original parameter values within ParmEstAbsTol/ParmEstRelTol.');
            
            % Reinstate the original parameter values:
            testCase.Dist.ResetParms(HoldParms);
        end % ProbitmAFCChiSqEstTest
        
    end  % methods(Test)
    
end

