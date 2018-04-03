function [Computed] = CupiTest(varargin)
% function CupiTest(dist,cto):  Exercise and check the distribution computations.
%
% Inputs:
% dist is an already-existing distribution or a string distribution name
% cto is a "CupiTest Options" control structure with fields listed & initialized in ctoDefaultLoad (so far):
%
% Outputs:
%   Computed: A structure with all of the stuff that is computed.
%   ctoDefault: A default control options structure.
%
% Special cases:
%   Call with no arguments to get a "CupiTest Options" control structure with the current default options.

% Copyright (C) 2018 Jeffrey Owen Miller
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program (00License.txt). If not, see 
%     <http://www.gnu.org/licenses/>.
%

%{
RunOptions
DistOptions
CombinedOptions = RunOptions, then override some based on DistOptions
%}

if nargin == 0
    Computed = ctoDefaultLoad;
    return
end

dist = varargin{1};

ctoDefault = ctoDefaultLoad;
if (nargin==1)||isempty(varargin{2})
    cto = ctoDefault;
else
    cto=varargin{2};
end

if cto.ConstructOnly
    fprintf(2,['CupiTest constructed distribution: ' dist.StringName '\n']);
    return
end

if cto.RandomParms
    dist.RandomParms;
end
fprintf(2,['CupiTest checks for distribution: ' dist.StringName '\n']);

if ~cto.SkipPlot
    PlotDens(dist);  drawnow; % Plot PDF & CDF
end

Computed.TotalError = 0;
Computed.NErrorChecks = 0;
Computed.ErrorState = zeros(cto.MaxErrorChecks,1);
Computed.ErrorName = cell(cto.MaxErrorChecks,1);
Computed.MinOfRange = InverseCDF(dist,0.15);
Computed.MaxOfRange = InverseCDF(dist,0.85);

if cto.MeanSDetc
    MeanSDetc(dist);
end
if cto.Moments
    Moments(dist)
end
if cto.PDFCDFHaz
    PDFCDFHaz(dist)
end
if cto.PDFIntvsCDF
    PDFIntvsCDF(dist)
end
if cto.MGFs
    MGFs(dist);
end
if cto.FullIntegrals
    [Computed.RawInt, Computed.CenInt] = Integrals(dist,'Full Range',dist.LowerBound,dist.UpperBound);
end
if cto.PartialIntegrals
    [Computed.RawIntPartial, Computed.CenIntPartial] = Integrals(dist,'Partial Range',Computed.MinOfRange,Computed.MaxOfRange);
end

ConsistencyChecks;

cto.RandVals = Random(dist,cto.NRands,1);
if cto.CheckRands
    CheckRands(dist);
end

if (dist.NDistParms > 0) && (cto.AnyEst)
    ParameterEstimates(dist);
end
ErrorStates = Computed.ErrorState(1:Computed.NErrorChecks);
ErrorNames = Computed.ErrorName(1:Computed.NErrorChecks);
Computed.TotalError = sum(ErrorStates);
if Computed.TotalError > 0
    fprintf('%i errors found processing %s.\n',Computed.TotalError,dist.StringName);
    %    FailedErrors = '';
    errorpos = find(ErrorStates>0);
    s = ErrorNames(errorpos);
    FailedErrors = reshape(s,numel(s),1);
    fprintf('%s\n',FailedErrors{:});
else
    FailedErrors = ['CupiTest found no errors processing ' dist.StringName '.  :)'];
end

    function cto=ctoDefaultLoad()
        % These fields control program flow, e.g. which tests to perform
        cto.ConstructOnly = false;
        cto.SkipPlot = false;
        cto.SavePlot = 0;  % 0=never, 0.5=on error, 1=always
        cto.RandomParms = false;
        
        cto.MeanSDetc = true;
        cto.Moments = true;
        cto.PDFCDFHaz = true;
        cto.PDFIntvsCDF = true;
        cto.MGFs = true;
        % Chapter 5.7 of Numerical Recipes in C (http:%www.nrbook.com/a/bookcpdf/c5-7.pdf)
        % recommends a different (smaller) choice for cto.MGF_H, but I find their recommendation
        % is a bit too small for many distributions.
        cto.MGF_H = 1.0E-6; % Too large or too small produces numerical errors. This works quite well for normal.
        cto.FullIntegrals = true;
        cto.PartialIntegrals = true;
        cto.CheckRands = true;
        
        cto.AnyEst = true;
        cto.MLEEst = true;
        cto.MomentEst = true;
        cto.PercentileEst = true;
        cto.ChiSqEst = true;
        cto.ProbitYNMaxLikEst = true;
        cto.ProbitYNChiSqEst = true;
        cto.ProbitmAFCMaxLikEst = true;
        cto.ProbitmAFCChiSqEst = true;
        
        cto.XbyInverseErrorTolerance = 0.0001;
        cto.NPercentiles = 17;
        cto.Percentiles( 1) = 0.010;
        cto.Percentiles( 2) = 0.020;
        cto.Percentiles( 3) = 0.025;
        cto.Percentiles( 4) = 0.050;
        cto.Percentiles( 5) = 0.100;
        cto.Percentiles( 6) = 0.200;
        cto.Percentiles( 7) = 0.300;
        cto.Percentiles( 8) = 0.400;
        cto.Percentiles( 9) = 0.500;
        cto.Percentiles(10) = 0.600;
        cto.Percentiles(11) = 0.700;
        cto.Percentiles(12) = 0.800;
        cto.Percentiles(13) = 0.900;
        cto.Percentiles(14) = 0.950;
        cto.Percentiles(15) = 0.975;
        cto.Percentiles(16) = 0.980;
        cto.Percentiles(17) = 0.990;
        cto.NPowers = 5;
        cto.NThetas = 5;
        cto.ThetaStep = 0.02;
        cto.MatchToleranceAbs = 0.0001;
        cto.MatchToleranceRel = 0.001;
        cto.ErrorSignal = '****';
        cto.SmallErrorSignal = '**';
        cto.WantCDFIntegrals = true;
        cto.NRands = 1000;
        cto.NChiSqBins = 10;  % Used for Chi-square estimation & chi-square test of RNG.
        cto.NProbitBins = 10;
        cto.NProbitTrialsPerBin = 100;
        cto.ProbitmAFC = 2;  % The simplest example for the mAFC cases.
        cto.MaxErrorChecks = 200;
        
        % These options control printing:
        cto.PrintAll = true;
        %{
Possible Further options:
    Save PlotDens (only if problem?)
    Set for NONrandom starts
    newonly         Abort if the output file already exists.');
    redoonly        Abort if the output file does NOT already exist.');
    overwrite       Automatically overwrite existing output file.');
    nrands N        N of random numbers to generate.');
    nbins N         N of bins to use.');
    bintype         0, 1, 2, or 3.');
    nrngsamples N   N of random samples for testing RNG.');
    nrngbins N      N of data bins for testing RNG.');
    pctilecode N    Code for percentile values to generate:');
             1 ->   .05, .95');
             2 ->   .01, .05, .95, .99');
             3 ->   .01, .02, .025, .05, .95, .975, .98, .99');
             4 ->   .01, .02, .025, .05, .1, .2, .. .9, .95, .975, .98, .99');
    searchminstep m Minimum step size at which to terminate parm search.');
    targetcdf P     Target CDF value for percentile-based estimation.');

  Precision options should be set-able in distribution constructors
    integralprecision eps     Precision value for convergence of numerical integration.');
    inverseprecisionp eps     Precision value, in P, for computation of InverseCDF.');
    inverseprecisionx eps     Precision value, in X, for computation of InverseCDF.');
        %}
    end

    function MeanSDetc(dist)
        % This function just computes--no checks here.
        StringOut('Table of Summary Parameters:');
        Computed.TrueMean = Mean(dist);
        NumOut('Mean',Computed.TrueMean);
        Computed.TrueMedian = Median(dist);
        NumOut('Median',Computed.TrueMedian);
        Computed.TrueVariance = Variance(dist);
        NumOut('Variance',Computed.TrueVariance);
        Computed.TrueSD = SD(dist);
        NumOut('SD',Computed.TrueSD);
        if (Computed.TrueMean ~= 0)
            Computed.TrueCV = CV(dist);
            NumOut('CV',Computed.TrueCV);
        end
        Computed.TrueRawSkewness = RawSkewness(dist);
        NumOut('RawSkewness',Computed.TrueRawSkewness);
        Computed.TrueRelSkewness = RelSkewness(dist);
        NumOut('RelSkewness',Computed.TrueRelSkewness);
        Computed.TrueKurtosis = Kurtosis(dist);
        NumOut('Kurtosis',Computed.TrueKurtosis);
        Computed.TrueMinimum = Minimum(dist);
        NumOut('Minimum',Computed.TrueMinimum);
        Computed.TrueMaximum = Maximum(dist);
        NumOut('Maximum',Computed.TrueMaximum);
        Computed.TrueMMMSD = MMMSD(dist);
        NumOut('MMMSD',Computed.TrueMMMSD);
        Computed.TruePctileSkew75 = PctileSkew(dist,0.75);
        NumOut('PctileSkew(.75)',Computed.TruePctileSkew75);
        Computed.TruePctileSkew90 = PctileSkew(dist,0.90);
        NumOut('PctileSkew(.90)',Computed.TruePctileSkew90);
        Computed.TruePctileSkew99 = PctileSkew(dist,0.99);
        NumOut('PctileSkew(.99)',Computed.TruePctileSkew99);
        StringOut('');
    end

    function Moments(dist)
        NumColWid = '%18';
        Num2ColWid = '%36';
        IntColWid = '%6';
        DP = '.4';
        StringOut('Table of Moments:');
        if cto.PrintAll
            fprintf([IntColWid 's' Num2ColWid 's' Num2ColWid 's\n'],' ','----Unconditional----','----Conditional----');
            fprintf([IntColWid 's' NumColWid 's' NumColWid 's' NumColWid 's' NumColWid 's' '\n'],'Moment','Raw','Central','Raw','Central');
        end
        FS1 = [NumColWid DP 'f'];
        for I = 0:cto.NPowers
            Ip1 = I + 1;
            cto.RawMom(Ip1) = RawMoment(dist,I);
            cto.CenMom(Ip1) = CenMoment(dist,I);
            CondRaw = ConditionalRawMoment(dist,Computed.MinOfRange,Computed.MaxOfRange,I);
            if (I == 0)
                CheckMatch('Total conditional probability vs 1',CondRaw,1);
            end
            CondCen = ConditionalCenMoment(dist,Computed.MinOfRange,Computed.MaxOfRange,I);
            if (I == 1)
                CheckMatch('Total conditional central integral with power 1 vs 0',CondCen,0);
            end
            if cto.PrintAll
                fprintf([IntColWid 'i' FS1  FS1 FS1 FS1 '\n'],I,cto.RawMom(Ip1),cto.CenMom(Ip1),CondRaw,CondCen);
            end
        end
        StringOut('');
    end

    function []=PDFCDFHaz(dist)
        PctColWidDP = '%6.3f';
        PctColWid = '%6s';
        ColWidDP = '%12.4f';
        ColWid = '%12s';
        XbyInverseCheck = 'Computing X by CDF then InverseCDF.';
        XbyInverseBigErrorTolerance = 10*cto.XbyInverseErrorTolerance;
        StringOut('Table of X Values at Percentiles and Functions of X:');
        S=sprintf([PctColWid ColWid ColWid ColWid ColWid ColWid],'Pctile','X','PDF','CDF','XbyInvCDF','Hazard');
        StringOut(S);
        ErrorSumSq = 0;
        SumX = 0;
        SumXX = 0;
        for I = 1:cto.NPercentiles
            ThisPct = cto.Percentiles(I);
            ThisX = InverseCDF(dist,ThisPct);
            SumX = SumX + ThisX;
            SumXX = SumXX + ThisX^2;
            ThisPDF = PDF(dist,ThisX);
            ThisCDF = CDF(dist,ThisX);
            ThisXbyInverse = InverseCDF(dist,ThisCDF);
            ThisHazard = Hazard(dist,ThisX);
            S=sprintf([PctColWidDP ColWidDP ColWidDP ColWidDP ColWidDP ColWidDP],...
                ThisPct,ThisX,ThisPDF,ThisCDF,ThisXbyInverse,ThisHazard);
            StringOut(S);
            ErrorSumSq = ErrorSumSq + (ThisX - ThisXbyInverse)^2;
        end
        MeanX = SumX / cto.NPercentiles;
        VarOfX = (SumXX - cto.NPercentiles * MeanX^2) / cto.NPercentiles;
        CheckMatch(XbyInverseCheck,ErrorSumSq,0);
        % Old version does different checks depending on distribution.
        %if VarOfX == 0
        %    CheckMatch(XbyInverseError,ErrorSumSq,0);
        %    end
        %elseif (ErrorSumSq/VarOfX > XbyInverseBigErrorTolerance)
        %    CheckMatch([cto.ErrorSignal XbyInverseError],ErrorSumSq/VarOfX,XbyInverseBigErrorTolerance);
        %elseif (ErrorSumSq/VarOfX > cto.XbyInverseErrorTolerance)
        %    CheckMatch([cto.SmallErrorSignal XbyInverseError],ErrorSumSq/VarOfX,XbyInverseErrorTolerance);
        %end
        StringOut('');
    end

    function []=PDFIntvsCDF(dist)
        if dist.DistType ~= 'c'
            StringOut('Distribution is not continuous, so skip checks of integrated PDF values versus CDFs.');
            return;
        end
        PctColWidDP = '%6.3f';
        PctColWid = '%6s';
        ColWid = '%12s';
        ColWidDP = '%12.4f';
        PDFIntegralError = ' Error: Integral of PDF does not equal difference of CDFs.';
        PDFIntTolerance = 0.0001;
        StringOut('Table of Integrated PDF Values Between Successive Percentiles:');
        StringOut(sprintf([PctColWid ColWid ColWid ColWid],'Pctile','X','PDFInt','dCDF'));
        ErrorFound = 0;
        LastX = dist.LowerBound;
        LastCDF = 0;
        if dist.DistType == 'c'
            SmallDelta = 0;
        else
            SmallDelta = 1.0e-6;
        end
        WorstError = 0;
        for I = 1:cto.NPercentiles
            ThisPct = cto.Percentiles(I);
            ThisX = InverseCDF(dist,ThisPct);
            ThisPDFInt = IntegralXToNxPDF(dist,LastX,ThisX,0);
            ThisCDF = CDF(dist,ThisX);
            ThisError = abs(ThisCDF-LastCDF-ThisPDFInt);
            if ThisError > WorstError
                WorstError = ThisError;
            end
            if (abs(ThisCDF - LastCDF - ThisPDFInt) > PDFIntTolerance)
                ErrorFound = 1;
                ThisErrSignal = [' ' cto.ErrorSignal];
            else
                ThisErrSignal = '';
            end
            StringOut(sprintf([PctColWidDP ColWidDP ColWidDP ColWidDP '%s'], ...
                ThisPct,ThisX,ThisPDFInt,ThisCDF-LastCDF,ThisErrSignal));
            LastX = ThisX + SmallDelta;
            LastCDF = ThisCDF;
        end
        CheckMatch('Integral of PDF vs CDF: WorstError',WorstError,0);
        if ErrorFound
            StringOut([cto.ErrorSignal PDFIntegralError]);
        end
        StringOut('');
    end

    function []=MGFs(dist)
        MGFColWidS = '%44s';
        MGFColWidDP = '%44.4f';
        if (cto.ThetaStep >= 0.001)
            ThetaDP = 3;
        else
            ThetaDP = 5;
        end
        ThetaColWid = ThetaDP + 3;
        ThetaColWidDP = ['%' num2str(ThetaColWid) '.' num2str(ThetaDP) 'f'];
        ThetaColWidS = ['%' num2str(ThetaColWid) 's'];
        StringOut('Table of MGFs:');
        StringOut(sprintf([ThetaColWidS MGFColWidS],'Theta','MGF(Theta)'));
        if (mod(cto.NThetas,2)==1)
            LowI = 1; % Test Theta's symmetrically around 0
        else
            LowI = 0;
        end
        if (cto.NThetas > 0)
            for I = LowI:cto.NThetas
                ThisTheta = (I - cto.NThetas/2 - LowI + 0.5) * cto.ThetaStep;
                ThisMGF = MGF(dist,ThisTheta);
                if (ThisTheta==0)
                    MGF0 = ThisMGF;
                end
                StringOut(sprintf([ThetaColWidDP MGFColWidDP],ThisTheta,ThisMGF));
            end
        end
        MGFPlusH = MGF(dist,cto.MGF_H);
        MGFMinusH = MGF(dist,-cto.MGF_H);
        cto.MeanByMGF = (MGFPlusH - MGFMinusH) / (2*cto.MGF_H);
        cto.Mom2ByMGF = (MGFPlusH + MGFMinusH - 2*MGF0) / cto.MGF_H^2;
        StringOut('');
    end

    function [RawIntVals, CenIntVals] = Integrals(dist,Title,LowX, HighX)
        ColWid = 24;
        IntColWid = '%5';
        DP = 3;
        ColWids = ['%' num2str(ColWid) 's'];
        ColWidDPf = ['%' num2str(ColWid) '.' num2str(DP) 'f'];
        StringOut(['Table of ' Title ' Integrals:']);
        StringOut(sprintf([IntColWid 's' ColWids ColWids ColWids],'Power','Raw','Central','CDF'));
        RawIntVals = zeros(1,cto.NPowers+1);
        CenIntVals = zeros(1,cto.NPowers+1);
        for I = 0:cto.NPowers
            Ip1 = I+1;
            RawIntVals(Ip1) = IntegralXToNxPDF(dist,LowX,HighX,I);
            CenIntVals(Ip1) = IntegralX_CToNxPDF(dist,LowX,HighX,Computed.TrueMean,I);
            s = sprintf([IntColWid 'i ' ColWidDPf ' ' ColWidDPf],I,RawIntVals(Ip1),CenIntVals(Ip1));
            if cto.WantCDFIntegrals
                IntCDF = IntegralCDF(dist,LowX,HighX,I);
                s = [s sprintf(ColWidDPf,IntCDF)];
            end
            StringOut(s);
        end
        StringOut('');
    end

    function ConsistencyChecks()
        cto.ConsistencyOK = true;
        
        % Array indexing starts from 1 for 0th moment, etc.
        if cto.Moments;	CheckMatch('RawMoment(0)=1',cto.RawMom(1),1);	end;
        if cto.Moments;	CheckMatch('CenMoment(0)=1',cto.CenMom(1),1);	end;
        if cto.Moments&&cto.MeanSDetc;	CheckMatch('RawMoment(1)=Mu',cto.RawMom(2),Computed.TrueMean);	end;
        if cto.Moments;	CheckMatch('CenMoment(1)=0',cto.CenMom(2),0);	end;
        
        if cto.FullIntegrals;	CheckMatch('RawIntegral(0)=1',Computed.RawInt(1),1);	end;
        if cto.FullIntegrals;	CheckMatch('CenIntegral(0)=1',Computed.CenInt(1),1);	end;
        
        if cto.FullIntegrals&&cto.MeanSDetc;	CheckMatch('RawIntegral(1)=Mu',Computed.RawInt(2),Computed.TrueMean);	end;
        if cto.FullIntegrals&&cto.MeanSDetc;	CheckMatch('1st deriv. MGF(0)=Mu',cto.MeanByMGF,Computed.TrueMean);	end;
        if cto.FullIntegrals;	CheckMatch('CenIntegral(1)=0',Computed.CenInt(2),0);	end;
        
        if cto.Moments&&cto.MeanSDetc;	CheckMatch('RawMoment(2)=Var+Mu^2',cto.RawMom(3),Computed.TrueVariance+Computed.TrueMean^2);	end;
        if cto.MGFs&&cto.MeanSDetc;	CheckMatch('RawIntegral(2)=Var+Mu^2',Computed.RawInt(3),Computed.TrueVariance+Computed.TrueMean^2); end;
        if cto.MGFs&&cto.MeanSDetc;	CheckMatch('2nd deriv. MGF(0)=Var+Mu^2',cto.Mom2ByMGF,Computed.TrueVariance+Computed.TrueMean^2);	end;
        if cto.Moments&&cto.MeanSDetc;	CheckMatch('CenMoment(2)=Var',cto.CenMom(3),Computed.TrueVariance);	end;
        if cto.FullIntegrals&&cto.MeanSDetc;	CheckMatch('CenIntegral(2)=Var',Computed.CenInt(3),Computed.TrueVariance);	end;
        
        if cto.FullIntegrals&&cto.Moments;	CheckMatch('RawMoment(3)=RawIntegral(3)',cto.RawMom(4),Computed.RawInt(4));	end;
        if cto.FullIntegrals&&cto.Moments;	CheckMatch('CenMoment(3)=CenIntegral(3)',cto.CenMom(4),Computed.CenInt(4));	end;
        
        if cto.FullIntegrals&&cto.Moments;	CheckMatch('RawMoment(4)=RawIntegral(4)',cto.RawMom(5),Computed.RawInt(5));	end;
        if cto.FullIntegrals&&cto.Moments;	CheckMatch('CenMoment(4)=CenIntegral(4)',cto.CenMom(5),Computed.CenInt(5));	end;
        
        if cto.MeanSDetc;	CheckMatch('Rel vs Abs Skewness',Computed.TrueRawSkewness^3,Computed.TrueRelSkewness*Computed.TrueVariance*Computed.TrueSD);	end;
        if cto.Moments&&cto.MeanSDetc
            if cto.CenMom(4) > 0
                Temp = cto.CenMom(4)^(1/3);
            elseif cto.CenMom(4) < 0
                Temp = -(-cto.CenMom(4)^(1/3));
            else
                Temp = 0;
            end
            CheckMatch('Abs Skewness vs CenMom(3)',Computed.TrueRawSkewness,Temp);
            if (cto.CenMom(3) == 0)
                Temp = 0;
            else
                Temp = cto.CenMom(5)/cto.CenMom(3)^2;
            end
            CheckMatch('Kurtosis vs CenMom(4)',Computed.TrueKurtosis,Temp);
        end
        
        if cto.ConsistencyOK
            StringOut('Passed moment-based consistency checks.');
        end
        StringOut('');
    end

    function CheckAssertion(CheckName, Assertion)
        Computed.NErrorChecks = Computed.NErrorChecks + 1;
        if Assertion
            Computed.ErrorName{Computed.NErrorChecks} = CheckName;
            Computed.ErrorState(Computed.NErrorChecks) = 0;
        else
            Computed.ErrorName{Computed.NErrorChecks} = [cto.ErrorSignal ' -- ' CheckName ': ASSERTION FAILED!'];
            Computed.ErrorState(Computed.NErrorChecks) = 1;
            disp(Computed.ErrorName{Computed.NErrorChecks});
        end
    end

    function CheckMatch(CheckName, Val1, Val2)
        Computed.NErrorChecks = Computed.NErrorChecks + 1;
        Computed.ErrorName{Computed.NErrorChecks} = CheckName;
        BigErrorFactor = 10.0;
        Dif = abs(Val1 - Val2);
        if (Dif <= cto.MatchToleranceAbs)
            return;
        end
        Max = abs(Val1);
        if (abs(Val2) > Max)
            Max = abs(Val2);
        end
        RelDif = Dif / Max;
        if (RelDif < cto.MatchToleranceRel)
            return;
        end
        cto.ConsistencyOK = false;
        RelDif = RelDif * 100;  % Convert to percentage
        ErrString = [cto.ErrorSignal ' WARNING: Failed consistency check ' CheckName  '.  '];
        if (RelDif < Dif)
            ErrString = [ErrString '% error = ' num2str(RelDif)];
        else
            ErrString = [ErrString 'Absolute error = ' num2str(Dif)];
        end
        ErrString = sprintf('%s\n%s   Compared values are: %5.3f %5.3f',ErrString,cto.SmallErrorSignal,Val1,Val2);
        
        disp(ErrString);
        Computed.ErrorName{Computed.NErrorChecks} = ErrString;
        Computed.ErrorState(Computed.NErrorChecks) = 1;
    end

    function CheckRands(dist)
        if numel(find(cto.RandVals<dist.LowerBound)) > 0
            StringOut([cto.ErrorSignal ' WARNING: Random values less than lowerbound of ' num2str(dist.LowerBound)]);
        end
        if numel(find(cto.RandVals>dist.UpperBound)) > 0
            StringOut([cto.ErrorSignal ' WARNING: Random values greater than upperbound of ' num2str(dist.UpperBound)]);
        end
        BinMax = MakeBinSet(dist,cto.NChiSqBins,false);
        BinProb = FindBinProbs(dist,BinMax);
        [obschisqval, obschisqp] = obschisq(cto.RandVals,BinMax,BinProb);
        StringOut(['RNG check yields ObsChiSq = ' num2str(obschisqval) ' with ' num2str(cto.NChiSqBins) ' bins and p = ' num2str(obschisqp)] );
        CheckAssertion('Random numbers passed chi-square bin test.',(~isnan(obschisqp))&&(obschisqp > .01))
        %if (isnan(obschisqp)) || (obschisqp < .01)
        %    StringOut([cto.ErrorSignal ' WARNING: RNG test resulted in very low p value.']);
        %else
        %    StringOut('Random numbers passed chi-square bin test.');
        %end
        StringOut('');
    end

    function ReportEstResult(dist, StartError, EndError)
        % Const TitleWid == 29;  ErrColWid == 30; ErrDP == 10;
        StringOut(['starting error = ' num2str(StartError)]);
        StringOut(['stopping error = ' num2str(EndError)]);
        if StartError < EndError
            StringOut([cto.ErrorSignal ' WARNING: ERROR WENT UP.']);
        end
        StringOut(['ending distribution = ' dist.StringName]);
    end

    function []=ParameterEstimates(dist)
        ThisParmCodes = dist.DefaultParmCodes;
        % Hold distribution parameters so that they can be restored to original values for start of each fit.
        HoldParms = ParmValues(dist);
        
        if cto.MLEEst
            StringOut('Example of MLE Parameter Estimation:');
            StartError = -LnLikelihood(dist,cto.RandVals);
            EstML(dist,cto.RandVals,ThisParmCodes);
            EndError = -LnLikelihood(dist,cto.RandVals);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('MLE Estimation Reduced Error',StartError>=EndError);
            % There is really no useful summary to print out here.
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        
        if cto.MomentEst
            StringOut('Example of Method-of-Moments Estimation:');
            ObsMoments = dist.MomentsFromScores(cto.RandVals);
            StartError = MomentError(dist,ObsMoments);
            EstMom(dist,ObsMoments,ThisParmCodes);
            EndError = MomentError(dist,ObsMoments);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('Method-of-Moments Estimation Reduced Error',StartError>=EndError);
            ColWid = 12;
            StringOut([sprintf('%-*s',2*ColWid,' ') sprintf('%*s',ColWid+1,'Observed') sprintf('%*s',ColWid+1,'Predicted')]);
            DP = 3;
            StringOut(sprintf('%-*s %*.*f %*.*f',2*ColWid,'Mean',ColWid,DP,ObsMoments(1),ColWid,DP,dist.Mean) );
            if dist.NDistParms > 1
                StringOut(sprintf('%-*s %*.*f %*.*f',2*ColWid,'Variance',ColWid,DP,ObsMoments(2),ColWid,DP,dist.Variance) );
            end
            if dist.NDistParms > 2
                StringOut(sprintf('%-*s %*.*f %*.*f',2*ColWid,'CenMoment(3)',ColWid,DP,ObsMoments(3),ColWid,DP,dist.CenMoment(3)) );
            end
            if dist.NDistParms > 3
                StringOut(sprintf('%-*s %*.*f %*.*f',2*ColWid,'CenMoment(4)',ColWid,DP,ObsMoments(4),ColWid,DP,dist.CenMoment(4)));
            end
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        if cto.PercentileEst
            StringOut('Example of Percentile-matching Estimation:');
            [SampleValues, TargetCDF, NPctiles] = dist.PercentilesFromScores(cto.RandVals)
            StartError = PercentileError(dist,SampleValues,TargetCDF);
            EstPctile(dist,SampleValues,TargetCDF,ThisParmCodes);
            EndError = PercentileError(dist,SampleValues,TargetCDF);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('Percentile-matching Estimation Reduced Error',StartError>=EndError);
            ColWid = 20;
            DP = 10;
            StringOut(sprintf('%*s %*s %*s',ColWid,'X value',ColWid,'Observed CDF',ColWid,'Predicted CDF'));
            for I = 1:NPctiles
                StringOut(sprintf('%*.3f %*.*f %*.*f',ColWid,SampleValues(I),ColWid,DP,TargetCDF(I),ColWid,DP,CDF(dist,SampleValues(I))));
            end
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        if cto.ChiSqEst
            StringOut('Example of Chi-Square bin-based Estimation:');
            % Choose bins based on the original distribution.
            BinMax = MakeBinSet(dist,cto.NChiSqBins,false);
            BrainDeadHistc=histc(cto.RandVals,BinMax);
            BinProbs = BrainDeadHistc(1:numel(BinMax))/numel(cto.RandVals);
            StartError = GofFChiSq(dist,BinMax,BinProbs);
            EstChiSq(dist,BinMax,BinProbs,ThisParmCodes);
            EndError = GofFChiSq(dist,BinMax,BinProbs);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('Chi-Square bin-based Estimation Reduced Error',StartError>=EndError);
            ColWid = 14;
            StringOut(sprintf('%*s %*s %*s',2*ColWid,'BinTop',ColWid,'Observed Pr',ColWid,'Predicted Pr'));
            SumProb = 0;
            for IBin = 1:numel(BinProbs)
                SumProb = SumProb + BinProbs(IBin);
                if IBin>1
                    LowerProb = CDF(dist,BinMax(IBin-1));
                else
                    LowerProb = 0;
                end
                StringOut(sprintf('%*.*f %*.*f %*.*f',2*ColWid,DP,BinMax(IBin),...
                    ColWid,DP,BinProbs(IBin),ColWid,DP,CDF(dist,BinMax(IBin))-LowerProb) );
            end
            StringOut(sprintf('%*s %*.*f %*.*f',2*ColWid,['Above ' num2str(BinMax(end-1))],...
                ColWid,DP,1-SumProb,ColWid,DP,1-CDF(dist,BinMax(end-1))));
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        % Probit estimation starts here.
        
        GuessProb = 1 / cto.ProbitmAFC;
        NTrials = ones(cto.NProbitBins,1)*cto.NProbitTrialsPerBin;
        
        if cto.ProbitYNMaxLikEst || cto.ProbitYNMaxChiSq
            % Generate sample data for YN task:
            for iBin=1:cto.NProbitBins
                BinMax(iBin)=InverseCDF(dist,(iBin-0.5)/cto.NProbitBins);
                ThisP = CDF(dist,BinMax(iBin));
                NGreater(iBin) = binornd(NTrials(iBin),1-ThisP);
            end
        end
        
        if cto.ProbitYNMaxLikEst
            StringOut('Example of YNProbitLnLike Estimation:');
            StartError = -YNProbitLnLikelihood(dist,BinMax,NTrials,NGreater);
            EstProbitYNML(dist,BinMax,NTrials,NGreater,ThisParmCodes);
            EndError = -YNProbitLnLikelihood(dist,BinMax,NTrials,NGreater);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('YNProbitLnLike Estimation Reduced Error',StartError>=EndError);
            WriteObsPredYNProbit(dist,BinMax, NTrials, NGreater);
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        if cto.ProbitYNChiSqEst
            StringOut('Example of YNProbitChiSq Estimation:');
            StartError = YNProbitChiSq(dist,BinMax,NTrials,NGreater);
            EstProbitYNChiSq(dist,BinMax,NTrials,NGreater,ThisParmCodes);
            EndError = YNProbitChiSq(dist,BinMax,NTrials,NGreater);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('YNProbitChiSq Estimation Reduced Error',StartError>=EndError);
            WriteObsPredYNProbit(dist,BinMax, NTrials, NGreater);
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        if cto.ProbitmAFCMaxLikEst || cto.ProbitmAFCMaxChiSq
            % Generate sample data for cto.ProbitmAFC task:
            for iBin=1:cto.NProbitBins
                BinMax(iBin)=InverseCDF(dist,(iBin-0.5)/cto.NProbitBins);
                ThisP = CDF(dist,BinMax(iBin));
                ThisP = GuessProb + (1 - GuessProb) * ThisP;
                NGreater(iBin) = binornd(NTrials(iBin),1-ThisP);
            end
        end
        
        if cto.ProbitmAFCMaxLikEst
            StringOut('Example of mAFCProbitLnLike Estimation:');
            StartError = -mAFCProbitLnLikelihood(dist,cto.ProbitmAFC,BinMax,NTrials,NGreater);
            EstProbitmAFCML(dist,cto.ProbitmAFC,BinMax,NTrials,NGreater,ThisParmCodes);
            EndError = -mAFCProbitLnLikelihood(dist,cto.ProbitmAFC,BinMax,NTrials,NGreater);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('mAFCProbitLnLike Estimation Reduced Error',StartError>=EndError);
            WriteObsPredAFCProbit(dist,BinMax, NTrials, NGreater);
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
        if cto.ProbitmAFCChiSqEst
            StringOut('Example of mAFCProbitChiSq Estimation:');
            StartError = mAFCProbitChiSq(dist,cto.ProbitmAFC,BinMax,NTrials,NGreater);
            EstProbitmAFCChiSq(dist,cto.ProbitmAFC,BinMax,NTrials,NGreater,ThisParmCodes);
            EndError = mAFCProbitChiSq(dist,cto.ProbitmAFC,BinMax,NTrials,NGreater);
            ReportEstResult(dist, StartError, EndError);
            CheckAssertion('mAFCProbitChiSq Estimation Reduced Error',StartError>=EndError);
            WriteObsPredAFCProbit(dist,BinMax, NTrials, NGreater);
            StringOut('');
            % Reinstate distribution parameters for start of each fit.
            ResetParms(dist,HoldParms);
        end
        
    end

    function StringOut(S)
        if cto.PrintAll
            fprintf('%s\n',S);
        end
    end

    function NumOut(Title, R)
        %ColWid = 12;
        %NameColWid = 16;
        %DP = 4;
        %(Title :NameColWid,' = ' num2str(R:ColWid:DP))
        if cto.PrintAll
            fprintf('%16s = %12.4f\n',Title,R);
        end
        %If cto.Debug Then Writeln(Output,Title:NameColWid,' == ',R:ColWid:DP);
    end


    function WriteObsPredYNProbit(dist,BinMax, NTrials, NGreater)
        ColWid = 14;
        DP = 2;
        StringOut(sprintf('%*s %*s %*s',2*ColWid,'Stimulus Value',ColWid,'Observed Pr',ColWid,'Predicted Pr'));
        for iBin = 1:cto.NProbitBins
            StringOut(sprintf('%*.*f %*.*f %*.*f',2*ColWid,DP,BinMax(iBin),...
                ColWid,DP,NGreater(iBin)/NTrials(iBin),ColWid,DP,1-CDF(dist,BinMax(iBin))));
        end
    end

    function WriteObsPredAFCProbit(dist,BinMax, NTrials, NGreater)
        ColWid = 14;
        DP = 2;
        GuessProb = 1 / cto.ProbitmAFC;
        StringOut(sprintf('%*s %*s %*s',2*ColWid,'Stimulus Value',ColWid,'Observed Pr',ColWid,'Predicted Pr'));
        for iBin = 1:cto.NProbitBins
            PredProb = CDF(dist,BinMax(iBin));
            PredProb = GuessProb + (1 - GuessProb) * PredProb;
            PredProb = 1 - PredProb;
            StringOut(sprintf('%*.*f %*.*f %*.*f',2*ColWid,DP,BinMax(iBin),...
                ColWid,DP,NGreater(iBin)/NTrials(iBin),ColWid,DP,PredProb));
        end
    end

end


