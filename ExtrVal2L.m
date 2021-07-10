classdef ExtrVal2L < dContinuous
    % Extreme Value Type II distribution given by Luce (1986, p 508) with parameters shape (his alpha) and scale>0 (his beta).
    % As shape gets larger, distribution becomes more normal.
    % Unlike Luce, here we divide by scale so that larger scale values give distributions with larger scores.
    
    properties(SetAccess = protected)
        shape, scale
    end
    
    methods
        
        function obj=ExtrVal2L(varargin)
            obj=obj@dContinuous('ExtrVal2L');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.SearchOptions.MaxFunEvals = 2000;
            obj.SearchOptions.MaxIter = 2000;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExtrVal2L:Constructor', ...
                        'ExtrVal2L constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.shape = newparmvalues(1);
            obj.scale = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newshape = ifelse(ParmCodes(1)=='f', obj.shape, 1.05*obj.shape);
            newscale = ifelse(ParmCodes(2)=='f', obj.scale, 0.95*obj.scale);
            obj.ResetParms([newshape newscale]);
        end
        
        function []=ReInit(obj)
            assert(obj.shape>0,'ExtrVal2L shape must be > 0.');
            assert(obj.scale>0,'ExtrVal2L scale must be > 0.');
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            XX = X(InBounds) / obj.scale;
            XXP = XX.^(-obj.shape);
            thispdf(InBounds) = obj.shape*XXP.*exp(-XXP)./(XX*obj.scale);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            XX = X(InBounds) / obj.scale;
            XX = XX.^(-obj.shape);
            thiscdf(InBounds) = exp(-XX);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, ~, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            XX = -log(P);
            XX = XX.^(-1/obj.shape);
            thisval = XX * obj.scale;
        end
        
    end  % methods
    
end  % class ExtrVal2L



