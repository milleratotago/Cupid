classdef TruncParent < dTransMono
    % Parent class for various types of truncation distributions.
    
    properties(SetAccess = protected)
        LowerCutoffX, UpperCutoffX,   % Distribution cutoffs in terms of scores
        LowerCutoffP, UpperCutoffP,   % Distribution cutoffs in terms of probs.
        UnconditionalP, XTolerance,
        IthValOffset
    end
    
    methods
        
        function obj=TruncParent(DistName,BasisDist)
            obj=obj@dTransMono(DistName,BasisDist);
            obj.PDFScaleFactorKnown = true;
            obj.TransReverses = false;
        end

        function []=NewCutoffs(obj,LowerX,UpperX)
            obj.LowerCutoffX = max(LowerX,obj.BasisRV.LowerBound);
            obj.UpperCutoffX = min(UpperX,obj.BasisRV.UpperBound);
            obj.LowerCutoffP = obj.BasisRV.CDF(LowerX);
            obj.UpperCutoffP = obj.BasisRV.CDF(UpperX);
            obj.UnconditionalP = obj.UpperCutoffP - obj.LowerCutoffP;
            if obj.BasisRV.DistType == 'd'
                obj.UnconditionalP = obj.UnconditionalP + obj.BasisRV.PDF(obj.LowerCutoffX);
            end
        end
            
        function MakeTables(obj)
            Included = (obj.BasisRV.DiscreteX>=obj.LowerCutoffX) & (obj.BasisRV.DiscreteX<=obj.UpperCutoffX);
            obj.DiscreteX = obj.BasisRV.DiscreteX(Included);
            obj.NValues = numel(obj.DiscreteX);
            obj.DiscreteXmin = obj.BasisRV.DiscreteXmin(Included);
            obj.DiscreteXmax = obj.BasisRV.DiscreteXmax(Included);
            obj.DiscretePDF = obj.BasisRV.DiscretePDF(Included) / obj.UnconditionalP;
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.DiscreteCDF(end) = 1;
            obj.SetBinEdges;
            obj.LowerBound = obj.DiscreteXmin(1);
            obj.UpperBound = obj.DiscreteXmax(end);
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX;
        end
        
        function thisval = PDFScaleFactor(obj,X)
            thisval = ones(size(X))/obj.UnconditionalP;
        end

        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.DistType=='d'
                thiscdf = CDF@dDiscrete(obj,X);
                return;
            end
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf = zeros(size(X));
            thiscdf(InBounds) = obj.BasisRV.CDF(X(InBounds));
            thiscdf(InBounds) = (thiscdf(InBounds) - obj.LowerCutoffP) / obj.UnconditionalP;
            thiscdf(X>obj.UpperBound) = 1;
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            switch obj.DistType
                case 'c'
                    PUntruncated = obj.LowerCutoffP + P * obj.UnconditionalP;
                    thisval = InverseCDF(obj.BasisRV,PUntruncated);
                case 'd'
                    thisval = InverseCDF@dDiscrete(obj,P);
            end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Random(obj.BasisRV,varargin{:});
            for i=1:numel(thisval)
                while (thisval(i) < obj.LowerCutoffX) || (thisval(i) > obj.UpperCutoffX)
                    thisval(i) = Random(obj.BasisRV);
                end
            end
        end
        
        function thisval=RawMoment(obj,I)
            thisval = ConditionalRawMoment(obj.BasisRV,obj.LowerBound,obj.UpperBound,I);
        end
        
    end  % methods
    
end  % class TruncParent



