classdef DblMon < dContinuous
    % Double Monomial distribution with parameters tzero>0, delta>0, epsilon>1
    % Reference: Luce (1986), page 510.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        tzero, delta, epsilon,   % The distribution's parameters should be listed here.
        Frac1, deltaPlusOne, EpsMinusOne, CDFAttzero, SecondTermAttzero
    end
    
    methods
        
        function obj=DblMon(varargin)   % Constructor
            obj=obj@dContinuous('DblMon');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.XNearlyZero = eps;
            obj.CDFNearlyZero = 1e-7;  % Heavy trimming due to long tails. e-4 serious numerical errors
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('DblMon:Constructor', ...
                        'DblMon constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.tzero = newparmvalues(1);
            obj.delta = newparmvalues(2);
            obj.epsilon = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newtzero   = ifelse(ParmCodes(1)=='f', obj.tzero,  1.01*obj.tzero);
            newdelta   = ifelse(ParmCodes(2)=='f', obj.delta,  1.01*obj.delta);
            newepsilon = ifelse(ParmCodes(3)=='f', obj.epsilon,1.01*obj.epsilon);
            obj.ResetParms([newtzero newdelta newepsilon]);
        end
        
        function []=ReInit(obj)
            assert(obj.tzero>0,'DblMon tzero must be > 0.');
            assert(obj.delta>0,'DblMon delta must be > 0.');
            assert(obj.epsilon>1,'DblMon epsilon must be > 1.');
            obj.deltaPlusOne = obj.delta + 1;
            obj.EpsMinusOne = obj.epsilon - 1;
            obj.Frac1 = obj.deltaPlusOne*(obj.EpsMinusOne) / ( (obj.delta+obj.epsilon)*obj.tzero );
            obj.LowerBound = obj.XNearlyZero;
            obj.UpperBound = 2*obj.tzero;
            obj.Initialized = true;
            obj.CDFAttzero = CDF(obj,obj.tzero);
            obj.SecondTermAttzero = SecondTerm(obj,obj.tzero);
            StillLooking = true;
            while StillLooking
                obj.UpperBound = obj.UpperBound * 2;
                obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
                StillLooking = CDF(obj,obj.UpperBound-1) < obj.CDFNearlyOne;
            end
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(0,Parms(1)) NumTrans.GT2Real(0,Parms(2)) NumTrans.GT2Real(1,Parms(3))] ;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(0,Reals(1)) NumTrans.Real2GT(0,Reals(2)) NumTrans.Real2GT(1,Reals(3))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            xFrac = X / obj.tzero;
            XX = zeros(size(X));
            BigX = (X>obj.tzero) & (X>=obj.LowerBound) & (X<=obj.UpperBound);
            SmallX = (X<=obj.tzero) & (X>=obj.LowerBound) & (X<=obj.UpperBound);
            XX(SmallX) = xFrac(SmallX).^obj.delta;
            XX(BigX) = xFrac(BigX).^(-obj.epsilon);
            %             xxsize = size(XX)
            %             pdfsize = size(thispdf(InBounds))
            thispdf(InBounds) = obj.Frac1 * XX(InBounds);
            %            for i=1:numel(X)
            %                if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %                    xFrac = X(i) / obj.tzero;
            %                    if X(i) <= obj.tzero
            %                        XX = xFrac^obj.delta;
            %                    else
            %                        XX = xFrac^(-obj.epsilon);
            %                    end
            %                    thispdf(i) = obj.Frac1 * XX;
            %                end
            %            end
        end
        
        function thisval=SecondTerm(obj,X)
            xFrac = X / obj.tzero;
            XX = xFrac.^(-obj.epsilon)/(-obj.epsilon+1);
            thisval = obj.Frac1 .* XX .* X;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for i=1:numel(X)
                if X(InBounds)
                    if X(i) <= obj.tzero
                        xFrac = X(i) / obj.tzero;
                        XX = xFrac^obj.delta/(obj.delta+1);
                        thiscdf(i) = obj.Frac1 * XX * X(i);
                    else
                        thiscdf(i) = obj.CDFAttzero + SecondTerm(obj,X(i)) - obj.SecondTermAttzero;
                    end
                end
            end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.epsilon <= 2
                %                 thisval = NaN; % Theoretically correct but not for truncated version implemented here.
                thisval = Mean@dContinuous(obj);
            else
                thisval = obj.tzero * (obj.deltaPlusOne) * (obj.EpsMinusOne) / ( (obj.delta + 2) * (obj.epsilon - 2) );
            end
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.epsilon <= 3
                %                 thisval = NaN; % Theoretically correct but not for truncated version implemented here.
                thisval = Variance@dContinuous(obj);
            else
                thisval = obj.tzero^2 * (obj.deltaPlusOne) * (obj.EpsMinusOne) / ( (obj.delta + 3) * (obj.epsilon - 3) ) - Mean(obj)^2;
            end
        end
        
        %         function thisval=Hazard(obj,x)  % Luce gives this it does not agree with the PDF/CDF values so I do not think it is correct.
        %             assert(obj.Initialized,UninitializedError(obj));
        %             thisval=zeros(size(x));
        %             for i=1:numel(x)
        %                 if x(i) > obj.tzero
        %                     thisval(i) = (obj.EpsMinusOne) / x(i);
        %                 else
        %                     xOvertzero = x(i) / obj.tzero;
        %                     Denom = 1/obj.EpsMinusOne + 1/obj.deltaPlusOne - xOvertzero^obj.deltaPlusOne/obj.deltaPlusOne;
        %                     thisval(i) = xOvertzero^obj.epsilon * 1/obj.tzero / Denom;
        %                 end
        %             end
        %         end
        
    end  % methods
    
end  % class DblMon

