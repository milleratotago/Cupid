classdef Laplace < dContinuous
    % Laplace(location,scale>0)
    
    properties(SetAccess = protected)
        location, scale,
        OneOverTwoScale
    end
    
    methods
        
        function obj=Laplace(varargin)   % Constructor
            obj=obj@dContinuous('Laplace');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Laplace:Constructor', ...
                        'Laplace constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.location = newparmvalues(1);
            obj.scale = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newloc   = ifelse(ParmCodes(1)=='f', obj.location, obj.location+1);
            newscale = ifelse(ParmCodes(2)=='f', obj.scale,   1.1*obj.scale);
            obj.ResetParms([newloc newscale]);
        end
        
        function []=ReInit(obj)
            assert(obj.scale>0,'Laplace scale must be > 0.');
            obj.Initialized = true;
            obj.OneOverTwoScale = 0.5 / obj.scale;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            MinAbsDev = - abs( (X(InBounds)-obj.location)/obj.scale );
            thispdf(InBounds) = obj.OneOverTwoScale * exp(MinAbsDev);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            Dev = (X-obj.location) / obj.scale;
            thiscdf(InBounds&(X<=obj.location)) = exp(Dev(InBounds&(X<=obj.location))) / 2;
            thiscdf(InBounds&(X> obj.location)) = 1 - exp(-Dev(InBounds&(X> obj.location)))/2;
            % for i=1:numel(X)
            %     Dev = (X(i)-obj.location) / obj.scale;
            %     if X(i) <= obj.location
            %         thiscdf(i) = exp(Dev) / 2;
            %     else
            %         thiscdf(i) = 1 - exp(-Dev)/2;
            %     end
            % end
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            Sign = 2 * ((P>.5)-0.5);
            OneMinusStuff = 1 - (2*P(InBounds)-1).*Sign;
            thisval(InBounds) = obj.location - obj.scale * log(OneMinusStuff) .* Sign;
            % for i=1:numel(P)
            %     if P(i) < obj.CDFNearlyZero
            %         P(i) = obj.CDFNearlyZero;
            %     elseif P(i) > obj.CDFNearlyOne
            %         P(i) = obj.CDFNearlyOne;
            %     end
            %     if P(i) > 0.5
            %         Sign = 1;
            %     else
            %         Sign = -1;
            %     end
            %     OneMinusStuff = 1 - (2*P(i)-1)*Sign;
            %     thisval(i) = obj.location - obj.scale * log(OneMinusStuff) * Sign;
            % end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.location;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 2 * obj.scale^2;
        end
        
        function thisval=RawSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0;
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 6;
        end
        
    end  % methods
    
end  % class Laplace
