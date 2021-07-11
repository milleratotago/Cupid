classdef Frechet < dContinuous
    % Frechet distribution with parameters shape>0, scale>0, minval
    % See https://en.wikipedia.org/wiki/Frechet_distribution
    % NEWJEFF:
    %   Still numerical errors in unit tests.
    %   Wikipedia also has Skewness (finite for shape>3) & Kurtosis (finite for shape>4)
    
    methods(Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)), NumTrans.GT2Real(eps,Parms(2)), Parms(3)];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)), NumTrans.Real2GT(eps,Reals(2)), Reals(3)];
        end
        
        function thispdf = frechpdf(x,shape,scale,minval)
            thiscdf = Frechet.frechcdf(x,shape,scale,minval);
            thispdf = shape/scale * ( (x-minval)/scale ).^(-1-shape) .* thiscdf;
        end
        
        function thiscdf = frechcdf(x,shape,scale,minval)
            thiscdf = exp( -((x-minval)/scale).^(-shape));
        end
        
        function thisicdf = frechicdf(p,shape,scale,minval)
            thisicdf = zeros(size(p));
            for iel=1:numel(p)
                f = @(x) Frechet.frechcdf(x,shape,scale,minval) - p(iel);
                thisicdf(iel) = fzero(f,[minval realmax]);
            end
        end
        
    end % methods(Static)
    
    properties(SetAccess = protected)
        shape, scale, minval
    end
    
    methods
        
        function obj=Frechet(varargin)
            obj=obj@dContinuous('Frechet');
            obj.NDistParms = 3;
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.CDFNearlyOne = 0.99999;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{1},varargin{2},varargin{3}]);
                otherwise
                    ME = MException('Frechet:Constructor', ...
                        'Frechet constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.shape = newparmvalues(1);
            obj.scale = newparmvalues(2);
            obj.minval = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newShape = ifelse(ParmCodes(1)=='f', obj.shape, 0.9*obj.shape);
            newScale = ifelse(ParmCodes(2)=='f', obj.scale, 0.9*obj.scale);
            newminval = ifelse(ParmCodes(3)=='f', obj.minval, 1.1*obj.minval);
            obj.ResetParms([newShape, newScale, newminval]);
        end
        
        function []=ReInit(obj)
            assert(obj.shape>0,'Frechet shape parameter must be > 0.');
            assert(obj.scale>0,'Frechet scale parameter must be > 0.');
            obj.LowerBound = obj.minval+eps;
            obj.UpperBound = Frechet.frechicdf(obj.CDFNearlyOne,obj.shape,obj.scale,obj.minval);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
%             thispdf(InBounds) = obj.shape / obj.scale * ((X(InBounds) - obj.minval) / obj.scale).^(-1-obj.shape) .* exp(  (-(X(InBounds) - obj.minval)/obj.scale).^-obj.shape  );
            thispdf(InBounds) = Frechet.frechpdf(X(InBounds),obj.shape,obj.scale,obj.minval);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
%             thiscdf = exp( -((X-obj.minval)/obj.scale).^(-obj.shape));
            thiscdf(InBounds) = Frechet.frechcdf(X(InBounds),obj.shape,obj.scale,obj.minval);
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.shape > 1
                thisval = obj.minval + obj.scale*gamma(1- 1/obj.shape);
            else
                thisval = inf;
            end
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.shape > 2
                thisval = obj.scale^2*( gamma(1- 2/obj.shape) - gamma(1- 1/obj.shape)^2 );
            else
                thisval = inf;
            end
        end
        
        function thisval=Median(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.minval + obj.scale * log(2)^(-1/obj.shape);
        end
        
    end  % methods
    
end  % class Frechet

