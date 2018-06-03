classdef ExpSum < dContinuous
    % ExpSum(rate1, rate2): Sum of two exponentials with different rates.
    
    properties(SetAccess = protected)
        rate1, rate2,
        StoredMult, StoredRatio
    end
    
    methods
        
        function obj=ExpSum(varargin)
            obj=obj@dContinuous('ExpSum');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExpSum:Constructor', ...
                        'ExpSum constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.rate1 = newparmvalues(1);
            obj.rate2 = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newrate1 = ifelse(ParmCodes(1)=='f', obj.rate1, 0.9*obj.rate1);
            newrate2 = ifelse(ParmCodes(2)=='f', obj.rate2, 1.2*obj.rate2);
            obj.ResetParms([newrate1 newrate2]);
        end
        
        function []=ReInit(obj)
            assert(obj.rate1>0,'ExpSum rate1 must be > 0.');
            assert(obj.rate2>0,'ExpSum rate2 must be > 0.');
            assert(obj.rate1~=obj.rate2,'ExpSum rates must be different.');
            obj.StoredMult = obj.rate1*obj.rate2 / (obj.rate2-obj.rate1);
            obj.StoredRatio = 1 / (obj.rate1 - obj.rate2);
            obj.LowerBound = obj.CDFNearlyZero;
            obj.UpperBound = (1/obj.rate1 + 1/obj.rate2) * 1000;
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
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thispdf(InBounds) = (exp(-obj.rate1*X(InBounds))-exp(-obj.rate2*X(InBounds))) * obj.StoredMult;
            % for i=1:numel(X)
            %     if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %         thispdf(i) = (exp(-obj.rate1*X(i))-exp(-obj.rate2*X(i))) * obj.StoredMult;
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(X>=obj.UpperBound) = 1;
            thiscdf(InBounds) = obj.StoredRatio * ...
                        ( obj.rate1 * (1 - exp(-obj.rate2*X(InBounds))) - obj.rate2 * (1 - exp(-obj.rate1*X(InBounds))) );
            % for i=1:numel(X)
            %     if X(i) <= obj.LowerBound
            %     elseif X(i) >= obj.UpperBound
            %         thiscdf(i) = 1;
            %     else
            %         thiscdf(i) = obj.StoredRatio * ...
            %             ( obj.rate1 * (1 - exp(-obj.rate2*X(i))) - obj.rate2 * (1 - exp(-obj.rate1*X(i))) );
            %     end
            % end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1 / obj.rate1 + 1 / obj.rate2;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1 / obj.rate1^2 + 1 / obj.rate2^2;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            RN1 = rand(varargin{:});
            RN2 = rand(varargin{:});
            thisval = -log(RN1) / obj.rate1 - log(RN2) / obj.rate2;
        end
        
    end  % methods
    
end  % class ExpSum


