classdef Triangular < dContinuous
    % Triangular(min,max) distribution, with peak at center.
    
    properties(SetAccess = protected)
        min, max,
        center, heightfac
    end
    
    methods
        
        function obj=Triangular(varargin)
            obj=obj@dContinuous('Triangular');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Triangular:Constructor', ...
                        'Triangular constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.min = newparmvalues(1);
            obj.max = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, make the bounds a bit wider.
            OldLower = obj.LowerBound;
            OldUpper = obj.UpperBound;
            BoundShift = 0.051 * (OldUpper - OldLower);
            NewLower = ifelse(ParmCodes(1)=='f',obj.LowerBound,OldLower-BoundShift);
            NewUpper = ifelse(ParmCodes(2)=='f',obj.UpperBound,OldUpper+BoundShift);
            obj.ResetParms([NewLower NewUpper]);
        end
        
        function []=ReInit(obj)
            assert(obj.min<obj.max,'Triangular parameters must be min<max.');
            obj.center = (obj.min + obj.max) / 2;
            obj.heightfac = 4 / (obj.max-obj.min)^2;
            obj.LowerBound = obj.min;
            obj.UpperBound = obj.max;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % It is tempting to enforce the restriction that Parms(2) has to be greater than Parms(1), like this:
            % Reals = [Parms(1) NumTrans.GT2Real(Parms(1),Parms(2))];
            % Unfortunately, this creates problems if you want to estimate Parms(1) while holding Parms(2) fixed.
            % It may be better to use a different parameterization, e.g. Parms(1) = center & Parms(2) = width,
            % so that these parameters are independent.  See TriangularCW.
            Reals = Parms;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = Reals;
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            LowerPart = (X>=obj.LowerBound) & (X<=obj.center);
            thispdf(LowerPart) = (X(LowerPart) - obj.min) * obj.heightfac;
            UpperPart = (X<=obj.UpperBound) & (X>obj.center);
            thispdf(UpperPart) = (obj.max - X(UpperPart)) * obj.heightfac;
            % for i=1:numel(X)
            %     if (X(i) <= obj.min) || (X(i) >= obj.max)
            %     elseif X(i) <= obj.center
            %         thispdf(i) = (X(i) - obj.min) * obj.heightfac;
            %     else
            %         thispdf(i) = (obj.max - X(i)) * obj.heightfac;
            %     end
            % end
        end% 
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            thiscdf(X>obj.UpperBound) = 1;
            LowerPart = (X>=obj.LowerBound) & (X<=obj.center);
            UpperPart = (X<=obj.UpperBound) & (X>obj.center);
            thiscdf(LowerPart) = 0.5*(X(LowerPart) - obj.min) .* obj.PDF(X(LowerPart));
            thiscdf(UpperPart) = 1 - 0.5*(obj.max - X(UpperPart)) .* obj.PDF(X(UpperPart));
            % for i=1:numel(X)
            %     if X(i) <= obj.min
            %     elseif X(i) >= obj.max
            %         thiscdf(i) = 1;
            %     else
            %         HerePDF = PDF(obj,X(i));
            %         if X(i) <= obj.center
            %             thiscdf(i) = 0.5 * (X(i) - obj.min) * HerePDF;
            %         else
            %             % X >= obj.center
            %             thiscdf(i) = 1 - 0.5 * (obj.max - X(i)) * HerePDF;
            %         end
            %     end
            % end
        end

        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            thisval = zeros(size(P));
            for i=1:numel(P)
                if P(i) <= 0.5
                    thisval(i) = sqrt(2 * P(i) / obj.heightfac) + obj.min;
                else  % P >= 0.5
                    thisval(i) = -sqrt(2 * (1-P(i)) / obj.heightfac) + obj.max;
                end
            end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.center;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (obj.max-obj.min)^2 / 24;
        end
        
        function thisval=RawSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0;
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 2.4;   % 12 / 5
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            Triangular1 = rand(varargin{:});
            Triangular2 = rand(varargin{:});
            thisval = (Triangular1 - Triangular2) * (obj.center - obj.min) + obj.center;
        end
        
    end  % methods
    
end  % class Triangular


