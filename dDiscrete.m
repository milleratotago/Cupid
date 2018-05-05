classdef dDiscrete < dGeneric  % NWJEFF: Not vectorized
    % dDiscrete is an abstract class for discrete distributions.
    % It provides an option for computing and storing a complete table of
    % all X's, PDF(X), and CDF(X) when the distribution is initialized.
    
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

    properties(SetAccess = protected)
        StoredTablesInitialized,
        StoredX, StoredPDF, StoredCDF
    end
    
    properties(SetAccess = public)
        UseStoredTables   % Set to true if you want to use stored X/PDF/CDF tables.
        Grain  % A graininess factor, multiplied by eps, when checking CDF equality
        % If this is false, X/PDF/CDF are computed "on the fly".
    end
    
    methods(Abstract)
        thisval=LegalValue(obj,X)
        thisval=NearestLegal(obj,X)
        thisval=nIthValue(obj,I)
        thisval=NextValue(obj,X)   % Return the next- larger discrete value after X
    end
    
    methods
        
        % All descendants of the dDiscrete class must override nIthValue.
        
        % In addition, all descendants must override either
        %    nPDF, nCDF,  or  MakeTables.
        
        % If a descendant is EVER used with computation on the fly
        % (i.e., with "UseStoredTables" is sometimes false),
        % that descendant must override either nPDF or nCDF.
        
        % Note that the generic version of MakeTables provided here calls nIthValue and nPDF.
        
        %       function thisval=nIthValue(obj,I)  % "native" IthValue of X computed directly ("on the fly").
        %           ME = MException([obj.ThisFamilyName ':nIthValue'], ...
        %               'function dDiscrete.nIthValue should have been overridden.');
        %           throw(ME);
        %       end
        
        function obj=dDiscrete(FamName)
            obj=obj@dGeneric(FamName);
            obj.DistType = 'd';
            obj.UseStoredTables = false;
            obj.StoredTablesInitialized = false;
            obj.Grain = 1;
        end
        
        function []=MakeTables(obj)
            obj.StoredX = zeros(obj.NValues,1);
            obj.StoredPDF = zeros(obj.NValues,1);
            obj.StoredCDF = zeros(obj.NValues,1);
            sump = 0;
            for ival = 1:obj.NValues
                thisx = nIthValue(obj,ival);
                obj.StoredX(ival) = thisx;
                obj.StoredPDF(ival) = nPDF(obj,thisx);
                sump = sump + obj.StoredPDF(ival);
                obj.StoredCDF(ival) = sump;
            end
            obj.StoredTablesInitialized = true;
        end
        
        function thispdf=nPDF(obj,X)       % "native" PDF of X computed directly ("on the fly").
            thispdf = CDF(obj,X) - CDF(obj,X-eps(X));
        end
        
        function thiscdf=nCDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf=zeros(size(X));
            for i=1:numel(X)
                if X(i)<obj.LowerBound
                elseif X(i)>=obj.UpperBound
                    thiscdf(i) = 1;
                else
                    thiscdf(i) = 0;
                    I = 0;
                    StillLooking = true;
                    while StillLooking
                        I=I+1;
                        LesserX = IthValue(obj,I);
                        if LesserX <= X(i)
                            thiscdf(i) = PDF(obj,LesserX) + thiscdf(i);
                        end
                        StillLooking = LesserX<X(i);
                    end
                end
            end
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.UseStoredTables
                if ~obj.StoredTablesInitialized
                    MakeTables(obj);
                end
                thispdf = zeros(size(X));
                for i=1:numel(X)
                    i_of_x = find(obj.StoredX==X(i),1);
                    if numel(i_of_x) == 0
                        thispdf(i) = 0;
                    else
                        thispdf(i) = obj.StoredPDF(i_of_x);
                    end
                end
            else
                thispdf = nPDF(obj,X);
            end
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.UseStoredTables
                if ~obj.StoredTablesInitialized
                    MakeTables(obj);
                end
                thiscdf = zeros(size(X));
                for i=1:numel(X)
                    i_of_x = find(obj.StoredX>X(i),1);
                    if numel(i_of_x) == 0
                        thiscdf(i) = 1;
                    elseif i_of_x == 1
                        thiscdf(i) = 0;
                    else
                        thiscdf(i) = obj.StoredCDF(i_of_x-1);
                    end
                end
            else
                thiscdf = nCDF(obj,X);
            end
        end
        
        function thisval=IthValue(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.UseStoredTables
                if ~obj.StoredTablesInitialized
                    MakeTables(obj);
                end
                thisval = zeros(size(I));
                for i=1:numel(I)
                    thisval(i) = obj.StoredX(I(i));
                end
            else
                thisval = nIthValue(obj,I);
            end
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            thisval=zeros(size(P));
            for iel=1:numel(thisval)
                Sum = 0;
                I = 0;
                StillLooking = true;
                while StillLooking
                    I = I+1;
                    X = IthValue(obj,I);
                    ThisPDF = PDF(obj,X);
                    Sum = ThisPDF + Sum;
                    StillLooking = (Sum+obj.Grain*eps(Sum) < P(iel)) && (I < obj.NValues);
                end
                thisval(iel) = X;
            end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval=rand(varargin{:});
            for iel=1:numel(thisval)
                sumpr=0;
                ix=0;
                StillLooking = true;
                while StillLooking
                    ix=ix+1;
                    v=IthValue(obj,ix);
                    sumpr=sumpr+PDF(obj,v);
                    StillLooking = (sumpr<thisval(iel)) & (ix<obj.NValues);
                end
                thisval(iel)=v;
            end
        end
        
        function [FromI, ToI]=IRange(obj,FromX,ToX)
            % Return the range of indices 1=smallest, etc, including
            % all values from FromX to ToX for a discrete random variable.
            assert(obj.Initialized,UninitializedError(obj));
            I = 0;
            StillLooking = true;
            while StillLooking
                I=I+1;
                StillLooking = IthValue(obj,I) < FromX;
            end
            FromI = I;
            I=I-1;
            StillLooking = true;
            while StillLooking
                I=I+1;
                StillLooking = (I <= obj.NValues) && (IthValue(obj,I) <= ToX);
            end
            ToI = I - 1;
        end
        
        function thisval=IntegralXToNxPDF(obj,FromX,ToX,N)
            % Returns the sum from FromX to ToX of X^N * PDF.   Note that the
            %  function value for N == 0 should be one and this property can
            %  be used as a check of the accuracy of the computation of PDF.
            thisval = 0;
            if FromX <= ToX
                [FromI, ToI]=IRange(obj,FromX,ToX);
                for I = FromI:ToI
                    X = IthValue(obj,I);
                    XtoN = X^N;
                    Pr = PDF(obj,X);
                    thisval = thisval + Pr * XtoN;
                end
            end
        end
        
        function thisval=IntegralX_CToNxPDF(obj,FromX,ToX,C,N)
            % Returns the sum from FromX to ToX of (X-C)^N * PDF
            % Note that the function value for N == 0 should be one and this property can
            % be used as a check of the accuracy of the computation of PDF.
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0;
            if FromX <= ToX
                [FromI, ToI]=IRange(obj,FromX,ToX);
                for I = FromI:ToI
                    X = IthValue(obj,I);
                    Pr = PDF(obj,X);
                    XtoN = (X-C)^N;
                    thisval = thisval + Pr * XtoN;
                end
            end
        end
        
        function thisval=ConditionalRawMoment(obj,FromX,ToX,I)
            assert(obj.Initialized,UninitializedError(obj));
            ConditionalP = CDF(obj,ToX) - CDF(obj,FromX) + PDF(obj,FromX);
            if (ConditionalP == 0)
                thisval = 0;
            else
                thisval = IntegralXToNxPDF(obj,FromX,ToX,I) / ConditionalP;
            end
        end
        
        function thisval=ConditionalCenMoment(obj,FromX,ToX,I)
            assert(obj.Initialized,UninitializedError(obj));
            ConditionalP = CDF(obj,ToX) - CDF(obj,FromX) + PDF(obj,FromX);
            if (ConditionalP == 0)
                thisval = 0;
            else
                ConditionalMu = ConditionalRawMoment(obj,FromX,ToX,1);
                thisval = IntegralX_CToNxPDF(obj,FromX,ToX,ConditionalMu,I) / ConditionalP;
            end
        end
        
        function thisval=IntegralCDF(obj,FromX,ToX,N)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0;
            for I = 1:obj.NValues
                ThisX = IthValue(obj,I);
                if (ThisX>=FromX) && (ThisX<=ToX)
                    thisval = thisval + ThisX^N * CDF(obj,ThisX);
                end
            end
        end
        
        function thisval=MGFrng(obj,Theta,FromX,ToX)
            assert(obj.Initialized,UninitializedError(obj));
            FromX = max(FromX,obj.LowerBound);
            ToX = min(ToX,obj.UpperBound);
            [FromI, ToI]=IRange(obj,FromX,ToX);
            thisval = 0;
            for I = FromI:ToI
                X = IthValue(obj,I);
                if (X >= FromX) && (X <= ToX)
                    thisval = thisval + exp(Theta*X) * PDF(obj,X);
                end
            end
        end
        
        function x = XsToPlot(obj)
            x = IthValue(obj,1:obj.NValues);
        end
        
        function BinMax=MakeBinSet(obj,MinPr)
            % This function creates an output row vector BinMax defining edges of bins covering the RV's range.
            % Each bin covers a probability of at least MinPr.  Edges are moved slightly above observations to
            % avoid numerical problems (e.g., with histcounts).
            PrUsed = 0;
            NSaved = 0;
            BinMax = zeros(1,obj.NValues);
            % ThisX = obj.LowerBound;
            while (PrUsed<1-MinPr)
                NSaved = NSaved + 1;
                ThisX = obj.InverseCDF(min(1,PrUsed+MinPr));
                PrUsed = obj.CDF(ThisX);
                BinMax(NSaved) = (9*ThisX + obj.NextValue(ThisX))/10;  % Move bin 10% of the distance to the next observation.
                % [NSaved, BinMax(NSaved), PrUsed, obj.UpperBound, obj.CDF(obj.UpperBound)]
            end
            BinMax(NSaved) = obj.UpperBound;
            BinMax(NSaved+1:end) = [];
        end
        
    end  % methods
    
end  % class dDiscrete


