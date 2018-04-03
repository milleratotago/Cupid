classdef ExpSumT < dContinuous
    % ExpSumT(rate1,rate2,cutoff) is the distribution of the sum of two exponentials with different rates,
    %  and the first exponential is truncated at a cutoff.
    % Rates must differ.
    
    properties(SetAccess = protected)
        rate1, rate2, cutoff,
        StoreDif, Store1, StoreMult
    end
    
    methods
        
        function obj=ExpSumT(varargin)
            obj=obj@dContinuous('ExpSumT');
            obj.NDistParms = 3;
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExpSumT:Constructor', ...
                        'ExpSumT constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.rate1 = newparmvalues(1);
            obj.rate2 = newparmvalues(2);
            obj.cutoff = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newrate1  = ifelse(ParmCodes(1)=='f', obj.rate1, 0.98*obj.rate1);
            newrate2  = ifelse(ParmCodes(2)=='f', obj.rate2, 1.02*obj.rate2);
            newcutoff = ifelse(ParmCodes(3)=='f', obj.cutoff, 0.99*obj.cutoff);
            obj.ResetParms([newrate1 newrate2 newcutoff]);
        end
        
        function []=ReInit(obj)
            assert(obj.rate1>0,'ExpSumT rate1 must be > 0.');
            assert(obj.rate2>0,'ExpSumT rate2 must be > 0.');
            assert(obj.cutoff>0,'ExpSumT cutoff must be > 0.');
            assert(obj.rate1~=obj.rate2,'ExpSumT rates must be different.');
            Exp1 = -(obj.rate1-obj.rate2)*obj.cutoff;
            Exp2 = -obj.rate1*obj.cutoff;
            obj.StoreDif = obj.rate1 - obj.rate2;
            obj.Store1 = 1 / (1 - exp(-obj.rate1*obj.cutoff));
            obj.StoreMult = obj.rate1 / (obj.rate1-obj.rate2) * (1 - exp(Exp1)) / (1 - exp(Exp2));
            obj.LowerBound = obj.XNearlyZero;
            obj.UpperBound = obj.cutoff + 1/obj.rate2 * 28; % 28 = Standard_Exponential.InverseCDF(CDFNearlyOne);
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))  NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(X>=obj.UpperBound) = 1;
            BelowCutoff = X <= obj.cutoff;
            thiscdf(InBounds & BelowCutoff) = obj.Store1 * (((1 - exp(-obj.rate1*X(InBounds & BelowCutoff)) ...
               - obj.rate1 / obj.StoreDif * exp(-obj.rate2*X(InBounds & BelowCutoff)) .* (( 1-exp(-obj.StoreDif*X(InBounds & BelowCutoff)) ))  )));
            thiscdf(InBounds & ~BelowCutoff) = 1 - exp(-obj.rate2*X(InBounds & ~BelowCutoff)) * obj.StoreMult;

            % for i=1:numel(X)
            %     if X(i) <= obj.LowerBound
            %     elseif X(i) >= obj.UpperBound
            %         thiscdf(i) = 1;
            %     elseif X(i) <= obj.cutoff
            %         thiscdf(i) =  obj.Store1 * (((1 - exp(-obj.rate1*X(i)) - obj.rate1 / obj.StoreDif * exp(-obj.rate2*X(i)) * (( 1-exp(-obj.StoreDif*X(i)) ))  )));
            %     else  % X(i) > obj.cutoff
            %         thiscdf(i) = 1 - exp(-obj.rate2*X(i)) * obj.StoreMult;
            %     end
            % end
        end

        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            rand1 = rand(varargin{:})*(1-exp(-obj.cutoff*obj.rate1));
            rand1trunc = -log(1-rand1)/obj.rate1;
            rand2 = -log(rand(varargin{:}))/obj.rate2;
            thisval = rand1trunc + rand2;
        end

    end  % methods
    
end  % class ExpSumT

