classdef dDiscrete < dGeneric  % NWJEFF: Not vectorized
    % dDiscrete is an abstract class for discrete distributions.
    % Every distribution is stored as a list of values with associated PDF and CDF.
    
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
        DiscreteX    % Sorted list of the X values in this distribution...
        DiscretePDF  % ... and their probabilities...
        DiscreteCDF  % ... and the cumulative probabilities <= X.
        % These values define the minimum and maximums for each bin.
        DiscreteXmin
        DiscreteXmax
        DiscreteCDFmax
        StoredTablesInitialized
    end
    
    properties(SetAccess = public)
        XGrain     % A "graininess" factor for X that multiplies times eps(X) in checking Y values against DiscreteX values.
        % Y is accepted as equal to DiscreteX(i) if it is within XGrain*eps(DiscreteX(i))
        CDFGrain   % A "graininess" factor for CDF.  P is accepted as equal to DiscreteCDF(i) if it is within CDFGrain*eps(DiscreteCDF(i))
        ExhaustiveBins  % Set to true if all X values should be considered as falling within a bin.
    end
    
    methods
        
        function obj=dDiscrete(FamName)
            obj=obj@dGeneric(FamName);
            obj.DistType = 'd';
            obj.StoredTablesInitialized = false;
            obj.XGrain = 2;
            obj.CDFGrain = 2;  % InverseCDF returns a too-small value where CDF is increasing slowly
            obj.ExhaustiveBins = false;
        end
        
        function []=ClearBeforeResetParmsD(obj)
            obj.Initialized = false;
            obj.StoredTablesInitialized = false;
        end
        
        function []=TrimTables(obj,MinPDF,MaxCDF)
            % Trim the tables to remove X values with too-small PDF or too-large CDF
            HighPos = find(obj.DiscreteCDF>=MaxCDF,1);
            if numel(HighPos)==0
                HighPos = numel(obj.DiscreteX);
            end
            try % NWJEFF
                Keep = (obj.DiscretePDF>=MinPDF) & (obj.DiscreteX<=obj.DiscreteX(HighPos));
            catch
                disp('Problem')
            end
            obj.DiscreteX = obj.DiscreteX(Keep);
            obj.DiscretePDF = obj.DiscretePDF(Keep);
            obj.DiscreteCDF = obj.DiscreteCDF(Keep);
            obj.DiscreteCDF(end) = 1;
            obj.NValues = numel(obj.DiscreteX);
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
        end
        
        function SetBinEdges(obj)
            % After the DiscreteX values have been determined, set the edges of the bins.
            % MakeTables should always call this at the end.
            if obj.ExhaustiveBins
                obj.DiscreteXmin = nan(size(obj.DiscreteX));
                obj.DiscreteXmax = nan(size(obj.DiscreteX));
                obj.DiscreteXmin(1) = -inf;
                obj.DiscreteXmax(end) = inf;
                obj.DiscreteXmin(2:end) = (obj.DiscreteX(1:end-1) + obj.DiscreteX(2:end)) / 2;
                obj.DiscreteXmax(1:end-1) = obj.DiscreteXmin(2:end);
            else
                thistol = obj.XGrain * eps(obj.DiscreteX);
                obj.DiscreteXmin = obj.DiscreteX - thistol;
                obj.DiscreteXmax = obj.DiscreteX + thistol;
            end
            obj.DiscreteCDFmax = min(1,obj.DiscreteCDF + obj.CDFGrain * eps(obj.DiscreteCDF));
            obj.StoredTablesInitialized = true;
        end
        
        function thisval=IthValue(obj,I)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisval = zeros(size(I));
            for i=1:numel(I)
                thisval(i) = obj.DiscreteX(I(i));
            end
        end
        
        function thisval=LegalValue(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisval = zeros(size(X));
            for i=1:numel(X)
                i_of_x = find(X(i)<=obj.DiscreteXmax,1);
                thisval(i) = X(i)>=obj.DiscreteXmin(i_of_x);
            end
        end
        
        function thisval=NearestLegal(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisval = zeros(size(X));
            for i=1:numel(X)
                absdiffs = abs(X(i)-obj.DiscreteX);
                [~, minpos] = min(absdiffs);
                thisval(i) = obj.DiscreteX(minpos);
            end
        end
        
        function thispdf=PDF(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thispdf = zeros(size(X));
            for i=1:numel(X)
                i_of_x = find(X(i)<=obj.DiscreteXmax,1);  % first bin such that X < DiscreteXmax
                if (numel(i_of_x) == 1) && ( X(i) >= obj.DiscreteXmin(i_of_x) )
                    thispdf(i) = obj.DiscretePDF(i_of_x);
                end
            end
        end
        
        function thiscdf=CDF(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thiscdf = zeros(size(X));
            for i=1:numel(X)
                i_of_x = find(X(i)<=obj.DiscreteXmax,1);
                if (numel(i_of_x) == 0)
                    thiscdf(i) = 1;
                elseif i_of_x > 1
                    thiscdf(i) = obj.DiscreteCDF(i_of_x);
                elseif X(i)>obj.DiscreteXmin(1)  % X is in the first bin
                    thiscdf(i) = obj.DiscreteCDF(1);
                else
                    thiscdf(i) = 0;      % X is below the first bin
                end
            end
        end
        
        function thisval=InverseCDF(obj,P)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            %           Useful for debugging:
            %           try
            %           if min(P)<0
            %               disp('NJeff too small');
            %           elseif max(P)>obj.DiscreteCDFmax(end)
            %               disp('NJeff too large');
            %           end
            %           catch
            %               disp('Problem here');
            %           end
            %           if ~(min(P)>=0)&&(max(P)<=obj.DiscreteCDFmax(end))
            %               disp('Combined fail');
            %           end
            if (min(P)<0) || (max(P)>obj.DiscreteCDFmax(end))
                error(['InverseCDF requires 0<=P<=1 but found min(P) and max(P) of ' num2str(min(P)) ' and ' num2str(max(P))]);
            end
            thisval=zeros(size(P));
            for iel=1:numel(thisval)
                i_of_p = find(obj.DiscreteCDFmax>=P(iel),1);
                thisval(iel) = obj.DiscreteX(i_of_p);
            end
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisrand = rand(varargin{:});
            thisval = zeros(varargin{:});
            for i=1:numel(thisrand)
                i_of_r = find(thisrand(i)<=obj.DiscreteCDF,1);
                thisval(i)=obj.DiscreteX(i_of_r);
            end
        end
        
        function [FromI,ToI]=IRange(obj,FromX,ToX)
            % Return the range of indices 1=smallest, etc, including
            % all values from FromX to ToX for a discrete random variable.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            FromI = find(FromX>=obj.DiscreteXmin,1,'last');
            if numel(FromI)==0
                FromI = 1;
            end
            ToI = find(ToX<=obj.DiscreteXmax,1);
            if numel(ToI)==0
                ToI = obj.NValues;
            end
        end
        
        function thisval=SumX_CToNxPDF(obj,FromI,ToI,C,N)
            thisval = 0;
            for i = FromI:ToI
                XtoN = (obj.DiscreteX(i)-C)^N;
                Pr = obj.DiscretePDF(i);
                thisval = thisval + Pr * XtoN;
            end
        end
        
        function thisval=IntegralXToNxPDF(obj,FromX,ToX,N)
            % Returns the sum from FromX to ToX of X^N * PDF.   Note that the
            %  function value for N == 0 should be one and this property can
            %  be used as a check of the accuracy of the computation of PDF.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisval = 0;
            if FromX <= ToX
                [FromI,ToI] = IRange(obj,FromX,ToX);
                thisval = obj.SumX_CToNxPDF(FromI,ToI,0,N);
            end
        end
        
        function thisval=IntegralX_CToNxPDF(obj,FromX,ToX,C,N)
            % Returns the sum from FromX to ToX of (X-C)^N * PDF
            % Note that the function value for N == 0 should be one and this property can
            % be used as a check of the accuracy of the computation of PDF.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisval = 0;
            if FromX <= ToX
                [FromI,ToI] = IRange(obj,FromX,ToX);
                thisval = obj.SumX_CToNxPDF(FromI,ToI,C,N);
            end
        end
        
        function thisval=ConditionalMoment(obj,FromX,ToX,Around,N)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.StoredTablesInitialized
                MakeTables(obj);
            end
            thisval = 0;
            if FromX > ToX
                return;
            end
            [FromI,ToI] = IRange(obj,FromX,ToX);
            ConditionalP = obj.DiscreteCDF(ToI) - obj.DiscreteCDF(FromI) + obj.DiscretePDF(FromI);
            if (ConditionalP == 0)
                thisval = 0;
            else
                thisval = obj.SumX_CToNxPDF(FromI,ToI,Around,N) / ConditionalP;
            end
        end
        
        function thisval=ConditionalRawMoment(obj,FromX,ToX,N)
            thisval = obj.ConditionalMoment(FromX,ToX,0,N);
        end
        
        function thisval=ConditionalCenMoment(obj,FromX,ToX,N)
            ConditionalMu = ConditionalMoment(obj,FromX,ToX,0,1);  % Conditional raw mean
            thisval = obj.ConditionalMoment(FromX,ToX,ConditionalMu,N);
        end
        
        function thisval=IntegralCDF(obj,FromX,ToX,N)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 0;
            if FromX > ToX
                return;
            end
            [FromI,ToI] = IRange(obj,FromX,ToX);
            for i = FromI:ToI
                XtoN = obj.DiscreteX(i)^N;
                Pr = obj.DiscreteCDF(i);
                thisval = thisval + Pr * XtoN;
            end
        end
        
        function thisval=MGFrng(obj,Theta,FromX,ToX)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 0;
            if FromX > ToX
                return;
            end
            % FromX = max(FromX,obj.LowerBound);
            % ToX = min(ToX,obj.UpperBound);
            [FromI, ToI]=IRange(obj,FromX,ToX);
            for i = FromI:ToI
                X = obj.DiscreteX(i);
                thisval = thisval + exp(Theta*X) * obj.DiscretePDF(i);
            end
        end
        
        function x = XsToPlot(obj)
            % x = IthValue(obj,1:obj.NValues);
            x = MakeBinSetD(obj,.01);
        end
        
        function [BinMax,BinProb]=MakeBinSetD(obj,MinPr)
            % This function creates an output row vector BinMax defining the top edges of bins covering the RV's range.
            % Each bin covers a probability of at least MinPr.  Edges are moved slightly above observations to
            % avoid numerical problems (e.g., with histcounts).
            if obj.NValues <= 20
                BinMax = obj.DiscreteX;
                BinProb = obj.DiscretePDF;
            else
                PrUsed = 0;
                NSaved = 0;
                
                % Find the indices of the X's at the top of each bin:
                BinIs = zeros(1,obj.NValues);  % Allocate the maximum space needed, if every value has its own bin.
                while (PrUsed<1-MinPr)
                    NSaved = NSaved + 1;
                    i_of_p = find(obj.DiscreteCDFmax>=PrUsed+MinPr,1);
                    if numel(i_of_p)==0
                        i_of_p = obj.NValues;  % Use the largest value
                    end
                    BinIs(NSaved) = i_of_p;
                    PrUsed = obj.DiscreteCDF(i_of_p);
                end
                if PrUsed<1
                    % NSaved = NSaved + 1;  % Add the remaining probability into the final bin
                    BinIs(NSaved) = numel(obj.DiscreteX);
                end
                BinIs(NSaved+1:end) = [];
                try
                    BinMax = ((obj.DiscreteX(BinIs) + obj.DiscreteXmax(BinIs))/2);  % Good place for bin boundary
                catch
                    disp('Error computing BinMax.');
                end
                BinCDF = obj.DiscreteCDF(BinIs);
                BinProb = diff([0 BinCDF]);
            end
        end
        
    end  % methods
    
end  % class dDiscrete


