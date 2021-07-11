classdef NegativeBinomial < dDiscrete
    % NegativeBinomial(N,P) distribution of the number of failures before the Nth success with success probability P.
    % Useful reference: Evans, Hastings, & Peacock (1993)  who call my N, their X, and my X, their Y
    
    properties(SetAccess = protected)
        N, P
        NM1, Q, PtoN
    end
    
    methods (Static)
        
       function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.Bounded2Real(0,1,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2Bounded(0,1,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=NegativeBinomial(varargin)
            obj=obj@dDiscrete('NegativeBinomial');
            obj.ParmTypes = 'ir';
            obj.DefaultParmCodes = 'fr';
            obj.NDistParms = 2;
            
            obj.PDFNearlyZero = 2*eps(0);
            
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('NegativeBinomial:Constructor', ...
                        'NegativeBinomial constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsD(obj);
            obj.N = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.P = newparmvalues(2);
            ReInit(obj);
        end
        
        function []=PerturbParms(obj,ParmCodes)
            obj.N = ifelse(ParmCodes(1)=='f',obj.N,obj.N+1);
            obj.P = ifelse(ParmCodes(2)=='f',obj.P,0.985*obj.P);  % Very sensitive to local minimum problems
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            assert((obj.P>0)&(obj.P<1),'NegativeBinomial P must be in range (0-1).');
            obj.NM1 = obj.N - 1;
            obj.PtoN = obj.P^obj.N;
            obj.Q = 1 - obj.P;
            
            obj.Initialized = true;
            MaybeUB = floor(obj.Mean) + 1;
            MaybeLB = MaybeUB - 1;
            SumP = 0;
            UBReached = false;
            LBReached = false;
            while (~UBReached || ~LBReached)
                if ~LBReached
                    LoP = obj.PMF(MaybeLB);
                end
                if ~UBReached
                    HiP = obj.PMF(MaybeUB);
                end
                LBReached = LoP == 0;
                UBReached = HiP == 0;
                SumP = SumP + LoP + HiP;
                if SumP >= obj.CDFNearlyOne
                    LBReached = true;
                    UBReached = true;
                else
                    if ~UBReached
                        MaybeUB = MaybeUB + 1;
                    end
                    if ~LBReached
                        if MaybeLB > 0
                            MaybeLB = MaybeLB - 1;
                        else
                            LBReached = true;
                            LoP = 0;
                        end
                    end
                end
            end
            
            obj.LowerBound = MaybeLB;
            obj.UpperBound =  MaybeUB;
            obj.Initialized = true;
            MakeTables(obj);
            TrimTables(obj,eps(0),1);
            SetBinEdges(obj);
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Pr = PMF(obj,k)  % Probability mass function for k failures
            Pr = nchoosek(obj.N+k-1,k) * obj.PtoN * obj.Q^k;
        end
        
         function thisval=LegalValue(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = zeros(size(X));
            for i=1:numel(X)
                if (abs(round(X(i))-X(i))<eps) && (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
                    thisval(i) = true;
                end
            end
        end
        
        function thisval=NearestLegal(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = zeros(size(X));
            for iel=1:numel(X)
                if X(iel) <= obj.LowerBound
                    thisval(iel) = obj.LowerBound;
                elseif X(iel) < obj.UpperBound
                    thisval(iel) = round(X);
                else
                    X(iel) = obj.UpperBound;
                end
            end
        end
        
        function []=MakeTables(obj)
            % Make tables of obj.DiscreteX, obj.DiscretePDF, & obj.DiscreteCDF for k failures
            k = obj.LowerBound:obj.UpperBound;
            obj.DiscreteX = k;
            obj.DiscretePDF = zeros(1,numel(k));
            obj.DiscreteCDF = zeros(1,numel(k));
            SumPr = 0;
            for i=1:numel(obj.DiscreteX)
                obj.DiscretePDF(i) = obj.PMF(k(i));
                SumPr = SumPr + obj.DiscretePDF(i);
                obj.DiscreteCDF(i) = SumPr;
            end
            obj.DiscreteCDF(end) = 1;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.N*obj.Q/obj.P;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.N*obj.Q/obj.P^2;
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (1 + obj.Q) / sqrt(obj.N*obj.Q);
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 3 + 6/obj.N + obj.P^2/(obj.N*obj.Q);
        end
        
        function thisval=MGF(obj,Theta)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.P^obj.N * (1 - obj.Q*exp(Theta))^(-obj.N);
        end
        
        %       function thisval=Random(obj,varargin)
        %       end
        
        %       function s=EstML(obj,Observations,varargin)
        %       end
        
    end  % methods
    
end  % class NegativeBinomial


