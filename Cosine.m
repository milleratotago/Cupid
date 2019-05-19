classdef Cosine < dContinuous
    % Cosine(mu,halfwidth)
    % See https://en.wikipedia.org/wiki/Raised_cosine_distribution
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu         % distribution mean
        halfwidth  % distribution extends from mu-halfwidth to mu+halfwidth
    end
    
    methods
        
        function obj=Cosine(varargin)
            obj=obj@dContinuous('Cosine');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.Initialized = true;
            switch nargin
                case 0
                    ResetParms(obj);
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Cosine:Constructor', ...
                        'Cosine constructor takes 0 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.halfwidth = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, move the mean 0.05*halfwidth farther from zero and increase halfwidth by 10%.
            MeanShift = 0.05 * obj.halfwidth;
            OldMean = obj.mu;
            if OldMean < 0
                NewMean = OldMean - MeanShift;
            else
                NewMean = OldMean + MeanShift;
            end
            NewMean = ifelse(ParmCodes(1)=='f', obj.mu, NewMean);
            NewHalfwidth = ifelse(ParmCodes(2)=='f', obj.halfwidth, 1.1*obj.halfwidth);
            obj.ResetParms([NewMean NewHalfwidth]);
        end
        
        function []=ReInit(obj)
            obj.LowerBound = obj.mu - obj.halfwidth;
            obj.UpperBound = obj.mu + obj.halfwidth;
            obj.Initialized = true;
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
            thispdf(InBounds) = 0.5/obj.halfwidth*(1+cos( (X(InBounds)-obj.mu)/obj.halfwidth * pi) );
            %            for i=1:numel(X)
            %                if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %                    thispdf(i) = 0.5/obj.halfwidth*(1+cos( (X(i)-obj.mu)/obj.halfwidth * pi) );
            %                end
            %             end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            z = (X(InBounds)-obj.mu)/obj.halfwidth;
            thiscdf(InBounds) = 0.5 * ( 1 + z + 1/pi*sin(z*pi) );
            %            for i=1:numel(X)
            %                if X(i) < obj.LowerBound
            %                elseif X(i) >= obj.UpperBound
            %                    thiscdf(i) = 1;
            %                else
            %                    z = (X(i)-obj.mu)/obj.halfwidth;
            %                    thiscdf(i) = 0.5 * ( 1 + z + 1/pi*sin(z*pi) );
            %                end
            %            end
        end
        
        function thisval=MGF(obj,Theta)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = pi^2*sinh(obj.halfwidth*Theta) / (obj.halfwidth*Theta*(pi^2+obj.halfwidth^2*Theta^2)) * exp(obj.mu*Theta);
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.mu;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.halfwidth^2 * (1/3 - 2/pi^2);
        end
        
        function thisval=RawSkewness(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval= 0;
        end
        
        function thisval=Kurtosis(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval= 3 + 6 * (90 - pi^4) / ( 5*(pi^2-6)^2 );
        end
        
    end  % methods
    
end  % class Cosine

