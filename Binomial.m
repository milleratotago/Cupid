classdef Binomial < dDiscrete   % NWJEFF: Not vectorized
    % Binomial(N,P) distribution with parameters N_of_trials, Pr_of_success.
    
    % Useful reference: Evans, Hastings, & Peacock (1993)
    % Ling 1992 discusses an F-based computation that may be superior to the Normal or Poisson approximations.
    % Note that BinomialP is derived from this. Functions added here must also be added (or over-ridden) there.
    
    properties(SetAccess = protected)
        N, P,
        Q, NP, NQ, NPQ,
        % The next variable stores the type of approximation (if any) in use.
        Approx,     % 0 = Compute directly / Do not approximate
        % 1 = Normal approx
        % 2 = Poisson approx for x
        % 3 = Poisson approx for N-x
        Standard_Normal, PoissonBasis
    end
    
    properties(SetAccess = public)
        BinomialExact,  % If true, always use exact computation--do not approximate
        % If BinomialExact is false, the next three cutoffs are used to decide when and how to approximate.
        LargeNCutoff, LargeCellCutoff
    end
    
    methods
        
        function obj=Binomial(varargin)
            obj=obj@dDiscrete('Binomial');
            obj.ParmTypes = 'ir';
            obj.DefaultParmCodes = 'fr';
            obj.NDistParms = 2;
            obj.UseStoredTables = false;
            obj.BinomialExact = false;  % false allows approximations in some cases.
            obj.LargeNCutoff = 50;      % approximate if N larger than this
            obj.LargeCellCutoff = 10;   % Use normal approximation if at least this many per cell.
            obj.Standard_Normal = Normal(0,1);
            SetZExtreme(obj.Standard_Normal,10);  % Limit binomial to within this many SDs of mean.
            obj.PoissonBasis = Poisson();  % Create the Poisson basis.
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Binomial:Constructor', ...
                        'Binomial constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.N = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.P = newparmvalues(2);
            ReInit(obj);
        end
        
        function []=PerturbParms(obj,ParmCodes)
            obj.N = ifelse(ParmCodes(1)=='f',obj.N,obj.N+2);
            obj.P = ifelse(ParmCodes(2)=='f',obj.P,0.95*obj.P);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            assert((obj.P>0)&(obj.P<1),'Binomial P must be in range (0-1).');
            obj.LowerBound = 0;
            obj.UpperBound =  obj.N;
            obj.Initialized = true;
            obj.NP = obj.N*obj.P;
            obj.Q = 1 - obj.P;
            obj.NQ = obj.N*obj.Q;
            obj.NPQ = obj.NP*obj.Q;
            obj.Approx = 0;
            obj.LowerBound = 0;
            obj.UpperBound = obj.N;
            if (~obj.BinomialExact) && (obj.N > obj.LargeNCutoff)
                if (obj.NP > obj.LargeCellCutoff) && (obj.NQ > obj.LargeCellCutoff)
                    obj.Approx = 1;   % NORMAL
                    HalfWidth = obj.Standard_Normal.ZExtreme * sqrt(obj.NPQ);
                    obj.LowerBound = max(0,floor(obj.NP-HalfWidth));
                    obj.UpperBound = min(obj.N,floor(obj.NP+HalfWidth)+1);
                elseif obj.P <= obj.Q
                    % NWJEFF: Check Poisson basis next 2 cases.
                    obj.Approx = 2;   % POISSON on X
                    obj.PoissonBasis.ResetParms(obj.N*obj.P);
                    obj.LowerBound = 0;
                    obj.UpperBound = obj.PoissonBasis.UpperBound;
                else
                    obj.Approx = 3;   % POISSON ON obj.N-X
                    obj.PoissonBasis.ResetParms(obj.N*obj.Q);
                    obj.LowerBound = obj.N - obj.PoissonBasis.UpperBound;
                    obj.UpperBound = obj.N;
                end
            end
            while obj.PDF(obj.LowerBound)<=eps(obj.PDFNearlyZero)
                obj.LowerBound = obj.LowerBound + 1;
            end
            while obj.PDF(obj.UpperBound)<=eps(obj.PDFNearlyZero)
                obj.UpperBound = obj.UpperBound - 1;
            end
            obj.NValues = round(obj.UpperBound - obj.LowerBound) + 1;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.Bounded2Real(0,1,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2Bounded(0,1,Reals(2))];
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
        
        function thisval=nIthValue(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(I)>0&&max(I)<=obj.NValues,'Requested value at nonexistent position')
            thisval = obj.LowerBound+I-1;
        end
        
        function []=MakeTables(obj)
            obj.StoredX = zeros(obj.NValues,1);
            obj.StoredPDF = obj.StoredX;
            obj.StoredCDF = obj.StoredX;
            for i=1:obj.NValues
                obj.StoredX(i) = obj.LowerBound - 1 + i;
            end
            switch obj.Approx
                case 0
                    SumPr = 0;
                    T1 = 1;
                    T2 = 1;
                    T3 = (1-obj.P)^(obj.N);
                    for i=1:obj.NValues
                        K = obj.StoredX(i);
                        %                       T1 = nchoosek(obj.N,K);
                        %                       T2 = obj.P^K;
                        %                       T3 = (1-obj.P)^(obj.N-K);
                        
                        obj.StoredPDF(i) = T1 * T2 * T3;
                        
                        % Prepare updated T's for the next iteration:
                        T1 = T1 * (obj.N - K) / (K + 1);
                        T2 = T2 * obj.P;
                        T3 = T3 / (1-obj.P);
                        
                        % Version from Carl Nettelblad.
                        % Converted from C++ retrieved from http://www.nettelblad.se/binomcalc.cpp
                        % on 2013-10-30 (and elaborated slightly after emails):
                        SumPr = SumPr + obj.StoredPDF(i);
                        obj.StoredCDF(i) = SumPr;
                    end
                case 1
                    PrevSum = 0;
                    for i=1:obj.NValues
                        obj.StoredCDF(i) = CDF(obj.Standard_Normal, (obj.StoredX(i)+0.5-obj.NP)/sqrt(obj.NPQ) );
                        obj.StoredPDF(i) = obj.StoredCDF(i) - PrevSum;
                        PrevSum = obj.StoredCDF(i);
                    end
                case 2
                    PrevSum = 0;
                    for i=1:obj.NValues
                        obj.StoredCDF(i) = CDF(obj.PoissonBasis,obj.StoredX(i));
                        obj.StoredPDF(i) = obj.StoredCDF(i) - PrevSum;
                        PrevSum = obj.StoredCDF(i);
                    end
                case 3
                    PrevSum = 0;
                    for i=1:obj.NValues
                        obj.StoredCDF(i) = CDF(obj.PoissonBasis,obj.N-1-obj.StoredX(i));
                        obj.StoredPDF(i) = obj.StoredCDF(i) - PrevSum;
                        PrevSum = obj.StoredCDF(i);
                    end
            end
        end
        
        function thisval=NextValue(obj,X)
            thisval = X + 1;
        end
        
        function thispdf=nPDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            switch obj.Approx
                case 0
                    thispdf = zeros(size(X));
                    for i=1:numel(X)
                        if LegalValue(obj,X(i))
                            if obj.P == 0
                                if X(i) == 0
                                    thispdf(i) = 1;
                                else
                                    thispdf(i) = 0;
                                end
                            elseif obj.P == 1
                                if X == obj.N
                                    thispdf(i) = 1;
                                else
                                    thispdf(i) = 0;
                                end
                            else
                                % T1 = nchoosek(obj.N,round(X(i)));
                                % T2 = exp(X(i)*log(obj.P));
                                % T3 = exp((obj.N-X(i))*log(1-obj.P));
                                % thispdf(i) = T1 * T2 * T3;
                                thispdf(i) = binopdf(X(i),obj.N,obj.P);  %  Statistics Toolbox
                            end
                        end
                    end
                case   {1, 2, 3}
                    T1 = nCDF(obj,X);
                    T2 = nCDF(obj,X-1);
                    thispdf = T1 - T2;
            end
        end
        
        function thiscdf=nCDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            X = floor(X);
            switch obj.Approx
                case 0
                    thiscdf = zeros(size(X));
                    for i=1:numel(X)
                        % PSum = 0;
                        % for IVal=0:X(i)
                        %     PSum=PSum+PDF(obj,IVal);
                        % end
                        % thiscdf(i) = PSum;
                        thiscdf(i) = binocdf(X(i),obj.N,obj.P);  %  Statistics Toolbox
                    end
                case 1
                    thiscdf = CDF(obj.Standard_Normal, (X+0.5-obj.NP)/sqrt(obj.NPQ) );
                case 2
                    thiscdf = CDF(obj.PoissonBasis,X);
                case 3
                    thiscdf = 1 - CDF(obj.PoissonBasis,obj.N-1-X);
            end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.NP;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.NPQ;
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (obj.Q - obj.P) / sqrt(obj.NPQ);
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 3 - 6/obj.N + 1/obj.NPQ;
        end
        
        function thisval=MGF(obj,Theta)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (obj.P * exp(Theta) + obj.Q)^obj.N;
        end

        function thisval=Random(obj,varargin)
            if numel(varargin)==0
                varargin{1} = 1;
            end
            unirands = rand(varargin{:},obj.N);
            thisvar01 = unirands<obj.P;
            thissize = size(thisvar01);
            thisval = squeeze(sum(thisvar01,numel(thissize)));
        end
        
        function s=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            meanObs = mean(Observations);
            ResetParms(obj,[obj.N meanObs/obj.N]);
            BuildMyName(obj);
            s=obj.StringName;
        end
        
    end  % methods
    
end  % class Binomial


