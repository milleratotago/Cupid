classdef TruncatedX < dTransOf1
    % TruncatedX(BasisRV,LowerCutoffX,UpperCutoffX) produces the truncated version of the BasisRV,
    % truncating between the two indicated X cutoffs.
    
    properties(SetAccess = protected)
        LowerCutoffX, UpperCutoffX,   % Distribution cutoffs in terms of scores
        LowerCutoffP, UpperCutoffP,   % Distribution cutoffs in terms of probs.
        UnconditionalP, XTolerance,
        IthValOffset
    end
    
    methods
        
        function obj=TruncatedX(varargin)
            obj=obj@dTransOf1('TruncatedX');
            obj.NTransParms = 2;
            obj.TransParmCodes = 'rr';
            switch nargin
                case 0
                case 3
                    BuildMyBasis(obj,varargin{1});
                    ResetParms(obj,[ParmValues(obj.BasisRV) varargin{end-1} varargin{end}]);
                otherwise
                    ME = MException('TruncatedX:Constructor', ...
                        'TruncatedX constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues(1:end-2));
            obj.LowerCutoffX = newparmvalues(end-1);
            obj.UpperCutoffX = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes);
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLower = ifelse(ParmCodes(end-1)=='f', obj.LowerCutoffX,obj.LowerCutoffX-0.1);
            NewUpper = ifelse(ParmCodes(end)=='f', obj.UpperCutoffX,obj.UpperCutoffX+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewLower NewUpper]);
        end
        
        function []=ReInit(obj)
            assert(obj.LowerCutoffX<obj.UpperCutoffX,'TruncatedX distribution requires lower cutoff X < upper cutoff X.');
            if obj.LowerCutoffX <= obj.BasisRV.LowerBound
                obj.LowerBound = obj.BasisRV.LowerBound;
            else
                obj.LowerBound = obj.LowerCutoffX;
            end
            if obj.UpperCutoffX >= obj.BasisRV.UpperBound
                obj.UpperBound = obj.BasisRV.UpperBound;
            else
                obj.UpperBound = obj.UpperCutoffX;
            end
            obj.NDistParms = obj.BasisRV.NDistParms + 2;
            obj.DefaultParmCodes = [obj.BasisRV.DefaultParmCodes 'ff'];
            obj.DistType = obj.BasisRV.DistType;
            obj.NValues = obj.BasisRV.NValues;
            obj.IthValOffset = 0;
            switch obj.DistType
                case 'd'
                    obj.XTolerance = obj.XNearlyZero;
                    while IthValue(obj.BasisRV,obj.NValues) > obj.UpperBound
                        obj.NValues = obj.NValues - 1;
                    end
                    StillLooking = true;
                    while StillLooking
                        obj.IthValOffset = obj.IthValOffset + 1;
                        StillLooking = IthValue(obj.BasisRV,obj.IthValOffset) < obj.LowerBound;
                    end
                    obj.IthValOffset = obj.IthValOffset - 1;
                    obj.NValues = obj.NValues - obj.IthValOffset;
                case 'c'
                    obj.XTolerance = 0;
            end
            obj.LowerCutoffP = CDF(obj.BasisRV,obj.LowerBound-obj.XTolerance);
            obj.UpperCutoffP = CDF(obj.BasisRV,obj.UpperBound);
            obj.UnconditionalP = obj.UpperCutoffP - obj.LowerCutoffP;
            assert(obj.UnconditionalP>0,'TruncatedX distribution must include probability > 0.');
            
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = ParmValues(obj,varargin)
            parmvals = [obj.BasisRV.ParmValues obj.LowerCutoffX obj.UpperCutoffX];
        end
        
        function x = XsToPlot(obj)
            switch obj.DistType
                case 'c'
                    x = XsToPlot@dContinuous(obj);
                case 'd'
                    x = XsToPlot@dDiscrete(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.LowerCutoffX obj.UpperCutoffX];
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end-1:end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end-1:end);
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX;
        end
        
        function thisval = PDFScaleFactor(obj,X)
            % This is never called because dTransOf1.PDF is overwritten,
            % but it must be defined to make the class concrete.
            thisval = 1;
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            for i=1:numel(X)
                if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
                    UT = PDF(obj.BasisRV,X(i));
                    thispdf(i) = UT / obj.UnconditionalP;
                end
            end
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            for i=1:numel(X)
                if X(i) < obj.LowerBound
                elseif X(i) >= obj.UpperBound
                    thiscdf(i) = 1;
                else
                    XP = CDF(obj.BasisRV,X(i));
                    thiscdf(i) = (XP - obj.LowerCutoffP) / obj.UnconditionalP;
                end
            end
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            PUntruncated = obj.LowerCutoffP + P * (obj.UpperCutoffP - obj.LowerCutoffP);
            thisval = InverseCDF(obj.BasisRV,PUntruncated);
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
        
        function thisval=NearestLegal(obj,X)
            thisval = NearestLegal(obj.BasisRV,X);
            if thisval < obj.LowerBound
                thisval = IthValue(obj.BasisRV,obj.IthValOffset+1);
            elseif thisval > obj.UpperBound
                thisval = IthValue(obj.BasisRV,obj.IthValOffset+obj.NValues);
            end
        end
        
        function thisval=LegalValue(obj,X)
            thisval = (X >= obj.LowerBound) && (X <= obj.UpperBound) && LegalValue(obj.BasisRV,X);
        end
        
        function thisval=nIthValue(obj,I)
            thisval = IthValue(obj.BasisRV,obj.IthValOffset+I);
        end
        
        function thisval=RawMoment(obj,I)
            thisval = ConditionalRawMoment(obj.BasisRV,obj.LowerBound,obj.UpperBound,I);
        end
        
    end  % methods
    
end  % class TruncatedX



