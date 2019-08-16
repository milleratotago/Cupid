classdef dContinuous < dGeneric   % Calls by reference

    % dContinuous: Generic continuous random variable.

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
    
    properties(Hidden)  % Properties seen only by class methods
    end
    
    properties(SetAccess = private)    % These properties can only be set by the methods of this class.
    end
    
    properties(SetAccess = public)
        CDFdelta        % Step size when computing PDF as the difference between the CDFs of two nearby X values.
        InverseCDFTol   % Tolerance passed to fzero by InverseCDF
        IntegrateOverP  % T/F.  If T, integrals work slower, defining steps by percentiles rather than X values.
        % This is more accurate for many long-tailed distributions.
        % Integral tolerances
        IntegralPDFAbsTol, IntegralPDFRelTol
        IntegralCDFAbsTol, IntegralCDFRelTol
        IntegralMGFAbsTol, IntegralMGFRelTol
        IntegralPDFXNAbsTol, IntegralPDFXNRelTol
        IntegralPDFXmuNAbsTol, IntegralPDFXmuNRelTol
        IntegralCDFXNAbsTol, IntegralCDFXNRelTol
    end
    
    methods(Abstract)
    end  % Abstract methods
    
    methods
        
        function obj=dContinuous(FamName)   % Constructor
            obj=obj@dGeneric(FamName);  % Inherited constructor
            obj.InverseCDFTol = 0.0001;  % Tolerance on how close the CDF value must be to the target.
            obj.CDFdelta = 0.0001;
            obj.DistType = 'c';
            obj.UseStoredCDFs = false;
            obj.HaveStoredCDFs = false;
            obj.DefaultNBinsOfX = 200;
            obj.IntegrateOverP = false;
            % Tolerances for numerical integration with integral command; initialize to MATLAB default values.
            % The integral is regarded as converged if at least one of the tolerances is satisfied.
            obj.IntegralPDFAbsTol = 1e-10;  % For integrating PDF
            obj.IntegralPDFRelTol = 1e-6;
            obj.IntegralCDFAbsTol = 1e-10;  % For integrating CDF
            obj.IntegralCDFRelTol = 1e-6;
            obj.IntegralMGFAbsTol = 1e-10;  % For integrating MGF
            obj.IntegralMGFRelTol = 1e-6;
            MaxN = 100;  % Assume we will never want a moment greater than the 100th!
            obj.IntegralPDFXNAbsTol = 1e-10 * ones(1,MaxN);     % For integrating PDF*X^N for N>=1
            obj.IntegralPDFXNRelTol = 1e-6  * ones(1,MaxN);
            obj.IntegralPDFXmuNAbsTol = 1e-10 * ones(1,MaxN);   % For integrating PDF*(X-mu)^N for N>=1
            obj.IntegralPDFXmuNRelTol = 1e-6  * ones(1,MaxN);
            obj.IntegralCDFXNAbsTol = 1e-10 * ones(1,MaxN);     % For integrating CDF*X^N for N>=1
            obj.IntegralCDFXNRelTol = 1e-6  * ones(1,MaxN);
        end
        
        function []=ClearBeforeResetParmsC(obj)
            obj.Initialized = false;
            obj.HaveSplinePDFs = false;
            obj.HaveSplineCDFs = false;
            obj.HaveSplineInvCDFs = false;
            obj.HaveStoredCDFs = false;
        end

        function thispdf=PDF(obj,X)
            % Compute the PDF as the numerical derivative of the CDF.
            % This is used to supply PDF for distributions where only the CDF is known analytically.
            % In the descendent distribution, define the PDF function as PDFfromCDF.
            % For such a distribution CDF must of course be defined!
            % Examples: Kolmogorov, Order
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            Upper = X(InBounds)+obj.CDFdelta/2;
            Upper(Upper>obj.UpperBound) = obj.UpperBound;
            Lower = X(InBounds)-obj.CDFdelta/2;
            Lower(Lower<obj.LowerBound) = obj.LowerBound;
            Dif = Upper - Lower;  % Not necessarily CDFdelta.
            thispdf(InBounds) = ( CDF(obj,Upper) - CDF(obj,Lower) ) ./ Dif;
        end
        
        function thiscdf=CDF(obj,X)
            % Get the CDF by numerical integral if it is not defined in the descendant.
            % This is an example of an anonymous function, integrating a function with
            % multiple parameters (obj & X), and integrating a class function.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thiscdf=zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(X>=obj.UpperBound) = 1;
            if obj.UseSplineCDF
                thiscdf(InBounds) = obj.GetSplineCDF(X(InBounds));
            else
                if obj.UseStoredCDFs && ~obj.HaveStoredCDFs
                    StoreCDFs(obj,obj.DefaultNBinsOfX);
                end
                if obj.UseStoredCDFs
                    for i=1:numel(X)
                        ilarger = find(obj.StoredXs>X(i),1);
                        if (numel(ilarger) == 0) || (ilarger == 1)
                            lowx = obj.LowerBound;
                            PreviousCDF = 0;
                        else
                            ismaller = ilarger - 1;
                            lowx = obj.StoredXs(ismaller);
                            PreviousCDF = obj.StoredCDFs(ismaller);
                        end
                        thiscdf(i)= PreviousCDF + integral(@(x) PDF(obj,x),lowx,X(i),'AbsTol',obj.IntegralPDFAbsTol,'RelTol',obj.IntegralPDFRelTol);
                    end
                else
                    for i=1:numel(X)
                        if X(i) < obj.LowerBound
                            thiscdf(i) = 0;
                        elseif X(i) > obj.UpperBound
                            thiscdf(i) = 1;
                        else
                            thiscdf(i)=integral(@(x) PDF(obj,x),obj.LowerBound,X(i),'AbsTol',obj.IntegralPDFAbsTol,'RelTol',obj.IntegralPDFRelTol);
                        end
                    end
                end
            end
        end
        
        function thisval=InverseCDF(obj,P)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if min(P)<0
                error(['InverseCDF requires all P>0 but this call has a P of ',num2str(min(P))]);
            end
            if max(P)>1
                error(['InverseCDF requires all P<1 but this call has a P of ',num2str(max(P))]);
            end
            thisval = zeros(size(P));
            for i=1:numel(P)
                try
                    % fzero fails if CDF is too big for LowerBound or too small for UpperBound
                    os=optimset('TolFun',obj.InverseCDFTol);
                    thisval(i) = fzero(@(x) CDF(obj,x)-P(i) , [obj.LowerBound obj.UpperBound], os);
                catch InverseCDFErr
                    if CDF(obj,obj.LowerBound)>P(i)
                        thisval(i) = obj.LowerBound;
                    else
                        thisval(i) = obj.UpperBound;
                    end
                end
            end
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = InverseCDF(obj,rand(varargin{:}));
        end
        
        function thisval=X_CToNxPDFfromP(obj,P,C,N)
            Xs = obj.InverseCDF(P);
            PDFs = obj.PDF(Xs);
            thisval =  PDFs .* (Xs-C).^N;
        end
        
        function thisval=EVFun(obj,Fun,FromX,ToX)
            % Expected value of any function of X.
            % This function returns the integral from FromX to ToX of any function Fun(x) * PDF(x).
            % It uses the distribution lower & upper bounds if FromX & ToX are not provided,
            % or if either is empty or outside the corresponding bound.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if nargin==2
                FromX = obj.LowerBound;
                ToX = obj.UpperBound;
            else
                if (numel(FromX)==0) || (FromX<obj.LowerBound)
                    FromX = obj.LowerBound;
                end
                if (numel(ToX)==0) || (ToX>obj.UpperBound)
                    ToX = obj.UpperBound;
                end
            end
            if FromX >= ToX
                thisval = 0;
            else
                thisval = integral(@(x) PDF(obj,x).*(Fun(x)), FromX,ToX); % ,'AbsTol',AbsTol,'RelTol',RelTol);
            end
        end
        
        function thisval=IntegralXToNxPDF(obj,FromX,ToX,N)
            % Returns the integral from FromX to ToX of X^N * PDF.   Note that the
            %  function value for N == 0 should be one and this property can
            %  be used as a check of the accuracy of the computation of PDF.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if N >= 1
                AbsTol = obj.IntegralPDFXNAbsTol(N);
                RelTol = obj.IntegralPDFXNRelTol(N);
            else
                AbsTol = obj.IntegralPDFAbsTol;
                RelTol = obj.IntegralPDFRelTol;
            end
            if FromX >= ToX
                thisval = 0;
            elseif obj.IntegrateOverP
                FromP = CDF(obj,FromX);
                ToP = CDF(obj,ToX);
                thisval=integral(@(p) X_CToNxPDFfromP(obj,p,0,N), FromP,ToP,'AbsTol',AbsTol,'RelTol',RelTol);
            else
                thisval=integral(@(x) PDF(obj,x).*(x.^N), FromX,ToX,'AbsTol',AbsTol,'RelTol',RelTol);
            end
        end
        
        function thisval=IntegralX_CToNxPDF(obj,FromX,ToX,C,N)
            % Returns the integral from FromX to ToX of (X-C)^N * PDF
            % Note that the function value for N == 0 should be one and this property can
            % be used as a check of the accuracy of the computation of PDF.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if N >= 1
                AbsTol = obj.IntegralPDFXmuNAbsTol(N);
                RelTol = obj.IntegralPDFXmuNRelTol(N);
            else
                AbsTol = obj.IntegralPDFAbsTol;
                RelTol = obj.IntegralPDFRelTol;
            end
            if FromX >= ToX
                thisval = 0;
            elseif obj.IntegrateOverP
                FromP = CDF(obj,FromX);
                ToP = CDF(obj,ToX);
                thisval=integral(@(p) X_CToNxPDFfromP(obj,p,C,N), FromP,ToP,'AbsTol',AbsTol,'RelTol',RelTol);
            else
                thisval=integral(@(x) PDF(obj,x).*((x-C).^N), FromX,ToX,'AbsTol',AbsTol,'RelTol',RelTol);
            end
        end
        
        function thisval=ConditionalRawMoment(obj,FromX,ToX,I)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            ConditionalP = CDF(obj,ToX) - CDF(obj,FromX);
            if (ConditionalP == 0)
                thisval = 0;
            else
                thisval = IntegralXToNxPDF(obj,FromX,ToX,I) / ConditionalP;
            end
        end
        
        function thisval=ConditionalCenMoment(obj,FromX,ToX,I)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            ConditionalP = CDF(obj,ToX) - CDF(obj,FromX);
            if (ConditionalP == 0)
                thisval = 0;
            else
                ConditionalMu = ConditionalRawMoment(obj,FromX,ToX,1);
                thisval = IntegralX_CToNxPDF(obj,FromX,ToX,ConditionalMu,I) / ConditionalP;
            end
        end
        
        function thisval=IntegralCDF(obj,FromX,ToX,N)
            % integrates CDF to find raw moments.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if N >= 1
                AbsTol = obj.IntegralCDFXNAbsTol(N);
                RelTol = obj.IntegralCDFXNRelTol(N);
            else
                AbsTol = obj.IntegralCDFAbsTol;
                RelTol = obj.IntegralCDFRelTol;
            end
            thisval=integral(@(x) CDF(obj,x).^N, FromX,ToX,'AbsTol',AbsTol,'RelTol',RelTol);
        end
        
        function thisval=MGFrng(obj,Theta,FromX,ToX)
            % Returns integral from FromX to ToX of exp(Theta*X) * PDF.
            % Note that MGF(obj, 0, LowerBound, UpperBound) should be one and this property
            %  can be used as a check of the accuracy of computing PDF.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            FromX = max(FromX,obj.LowerBound);
            ToX = min(ToX,obj.UpperBound);
            thisval=integral(@(x) exp(Theta*x).*PDF(obj,x),FromX,ToX,'AbsTol',obj.IntegralMGFAbsTol,'RelTol',obj.IntegralMGFRelTol);
        end
        
        function x = XsToPlot(obj)
            nsteps = 100;
            lowp = 0.001;
            highp = 0.999;
            lowx = InverseCDF(obj,lowp);
            highx = InverseCDF(obj,highp);
            stepsize = (highx-lowx)/nsteps;
            x = lowx:stepsize:highx;
        end
        
        function [BinMax,BinProb]=MakeBinSetC(obj,MinPr)
            % This function creates an output vector of length NBins defining NBins bins covering the RV's range.
            % The bottom of the first bin is implicitly obj.LowerBound, and the top of bin I is BinMax(I).
            % BinMax(NBins) always equals obj.UpperBound.
            % If EqualInX is true, the bins divide the range from LowerBound to UpperBound
            % into equal-length intervals.  If it is false, the bins are equal in probability.
            % This function only works for continuous distributions.
            NBins = ceil(1/MinPr);
            BinMax = zeros(1,NBins);
            BinMax(NBins) = obj.UpperBound;
            EqualInX = false;
            if EqualInX
                BinWidth = (obj.UpperBound - obj.LowerBound) / NBins;
                for iBin = 1:NBins-1
                    BinMax(iBin) = obj.LowerBound + BinWidth*(iBin);
                end
            else
                PropPerBin = 1 / NBins;
                for iBin = 1:NBins-1
                    BinMax(iBin) = InverseCDF(obj,PropPerBin*iBin);
                end
            end
            BinProb = obj.CDF(BinMax);
            BinProb = diff([0 BinProb]);
        end
        
        function Dmax = ksDmax(obj,x)
            % This fn is useful in computing the Kolmogorov-Smirnov test.
            % It computes the empirical CDF of the vector of observations in x, and it returns the
            % maximum difference between this empirical CDF and the theoretical CDF.
            sx = sort(x);
            n = numel(sx);
            Fnx = ( (1:n) / n )';  % Empirical Fn(x)
            Fx = obj.CDF(sx);
            dif = abs(Fnx - Fx);
            Dmax = max(dif);   % test statistic
        end
        
        function [p2tailed, Dmax] = kstest(obj,x,varargin)
            % Compute the Kolmogorov-Smirnov test to check whether the observations in x
            % are consistent with the current theoretical distribution in obj.
            % The optional varargin can be a KolmSmir distribution object.
            % If you will use the same KolmSmir distribution repeatedly then
            % you can speed up the code by pre-specifying it and passing it in.
            if numel(varargin) == 1
                ksdist = varargin{1};
            else
                N = numel(x);
                ksdist = KolmSmir(N);
            end
            Dmax = obj.ksDmax(x);
            p2tailed = 1 - ksdist.CDF(Dmax);
        end
        
    end  % methods
    
end  % class dContinRV

