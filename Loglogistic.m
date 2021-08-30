classdef Loglogistic < dContinuous
    % Loglogistic(scale>0,shape>0):  log(X) has a logistic distribution
    % Based on https://en.wikipedia.org/wiki/Log-logistic_distribution
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
      scale, shape  % called alpha & beta, respectively, in Wikipedia
      b  % s
    end
    
    properties(SetAccess = public)
        fsolveoptions  % Used by StartParmsMLE
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function parms = StartParmsMLE(X)
            % Median is scale, lower quartile
            pctiles = prctile(X,25:25:75);
            obs25pct = pctiles(1);
            obsmedian = pctiles(2);
            obs75pct = pctiles(3);
            estscale1 = -log(3)/log(obs25pct/obsmedian);
            estscale2 = log(3)/log(obs75pct/obsmedian);
            parms = [obsmedian, (estscale1+estscale2)/2];
            % Median formulas from Wikipedia
            % syms alpha beta obs25 obs75
            % solve(obs75 == 3^(1/beta) * alpha,beta)
            % solve(obs25 == 3^(-1/beta) * alpha,beta)
        end
        
    end
    
    methods
        
        function obj=Loglogistic(varargin)   % Constructor
            obj=obj@dContinuous('Loglogistic');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            obj.fsolveoptions = optimset('Display','off');
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Loglogistic:Constructor', ...
                        'Loglogistic constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.scale = newparmvalues(1);
            obj.shape = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            NewScale = ifelse(ParmCodes(1)=='f',obj.scale,1.05*obj.scale);
            NewShape = ifelse(ParmCodes(2)=='f',obj.shape,0.95*obj.shape);
            obj.ResetParms([NewScale NewShape]);
        end
        
        function []=ReInit(obj)
            if obj.scale<=0
                error('Loglogistic scale must be > 0.');
            end
            if obj.shape<=0
                error('Loglogistic shape must be > 0.');
            end
            obj.b = pi / obj.shape;
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = (obj.shape/obj.scale) * (X(InBounds)/obj.scale).^(obj.shape-1) ./ ...
                                ( 1 + (X(InBounds) / obj.scale).^obj.shape ).^2;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            Xpower = X(InBounds).^obj.shape;
            thiscdf(InBounds) = Xpower ./ (obj.scale^obj.shape + Xpower);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            thisval(InBounds) = obj.scale * (P(InBounds) ./ (1-P(InBounds))).^(1/obj.shape);
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.shape > 1
                thisval = (obj.scale*obj.b) / sin(obj.b);
            else
                warning('True mean undefined; approximating with mean of this truncated approximation distribution.');
                thisval = Mean@dContinuous(obj);
            end
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.shape > 2
                thisval = obj.scale^2 * (2*obj.b / sin(2*obj.b) - obj.b^2 / sin(obj.b)^2);
            else
                warning('True variance undefined; approximating with variance of this truncated approximation distribution.');
                thisval = Variance@dContinuous(obj);
            end
        end
        
%         function s = EstMom(obj,TargetVals)
%             % Based on the following moment functions from Wikipedia:
%             % (a = scale, b = shape)
%             % mu == a*b/sin(b)
%             % var == a^2 * (2*b / sin(2*b) - b^2 / sin(b)^2)
%             % Here is the array of functions, all of which are to be zero'ed
%             F = @(p) [...
%                 ( TargetVals(1) - p(1)*p(2)/sin(p(2))); ...
%                 ( TargetVals(2) - p(1)^2 * ( 2*p(2))/sin(2*p(2)) - p(2)^2 / sin(p(2))^2 ) ...
%                 ];
%             x0 = [10 2];
%             EndingVals = fsolve(F,x0,obj.fsolveoptions);
%             obj.ResetParms(EndingVals);
%             BuildMyName(obj);
%             s=obj.StringName;
%         end
        
    end  % methods
    
end  % class Loglogistic



