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
        ZExtreme  % Extreme Z value used in normal approximation
    end
    
    methods
        
        function obj=Binomial(varargin)
            obj=obj@dDiscrete('Binomial');
            obj.ParmTypes = 'ir';
            obj.DefaultParmCodes = 'fr';
            obj.NDistParms = 2;
            obj.BinomialExact = false;  % false allows approximations in some cases.
            obj.ZExtreme = 10;
            
            % Problems can arise if LowerBound and UpperBound are adjusted to exclude
            % some values with nonzero PDFs. In that case, transformations that reverse
            % the order of Xs (e.g., multiply by -1) may have incorrect CDFs due to
            % ignoring Xs outside the bounds.
            obj.PDFNearlyZero = 2*eps(0);
            
            obj.LargeNCutoff = 100;      % approximate if N larger than this
            obj.LargeCellCutoff = 20;   % Use normal approximation if at least this many per cell.
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
            assert((obj.P>0)&(obj.P<1),'Binomial P must be in range (0-1).');
            obj.LowerBound = 0;
            obj.UpperBound =  obj.N;
            obj.Initialized = true;
            obj.NP = obj.N*obj.P;
            obj.Q = 1 - obj.P;
            obj.NQ = obj.N*obj.Q;
            obj.NPQ = obj.NP*obj.Q;
            if obj.BinomialExact || (obj.N<=obj.LargeNCutoff)
                obj.Approx = 0;  % Exact
            elseif (obj.NP > obj.LargeCellCutoff) && (obj.NQ > obj.LargeCellCutoff)
                obj.Approx = 1;  % Normal
            elseif obj.P <= obj.Q
                obj.Approx = 2;   % POISSON on X
            else
                obj.Approx = 3;   % POISSON on N-X
            end
            MakeTables(obj);
            TrimTables(obj,eps(0),1);
            SetBinEdges(obj);
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
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
        
        function []=MakeTables(obj)
            switch obj.Approx
                case 0
                    [obj.DiscreteX, obj.DiscretePDF, obj.DiscreteCDF] = obj.MakeTableExact;
                case 1
                    [obj.DiscreteX, obj.DiscretePDF, obj.DiscreteCDF] = obj.MakeTableNormal;
                case {2, 3}
                    [obj.DiscreteX, obj.DiscretePDF, obj.DiscreteCDF] = obj.MakeTablePoisson;
            end
            obj.DiscreteCDF(end) = 1;
        end

        function [Xs, PDFs, CDFs] = MakeTableExact(obj)
            % Make tables of Xs, PDFs, & CDFs for 0-N successes
            SumPr = 0;
            T1 = 1;
            T2 = 1;
            T3 = (1-obj.P)^(obj.N);
            Xs = (0:obj.N);
            PDFs = zeros(1,obj.N+1);
            CDFs = zeros(1,obj.N+1);
            for i=1:numel(Xs)
                K = Xs(i);
                PDFs(i) = T1 * T2 * T3;
                % Prepare updated T's for the next iteration:
                T1 = T1 * (obj.N - K) / (K + 1);
                T2 = T2 * obj.P;
                T3 = T3 / (1-obj.P);
                % Version from Carl Nettelblad.
                % Converted from C++ retrieved from http://www.nettelblad.se/binomcalc.cpp
                % on 2013-10-30 (and elaborated slightly after emails):
                SumPr = SumPr + PDFs(i);
                CDFs(i) = SumPr;
            end
        end
        
        function [Xs, PDFs, CDFs] = MakeTableNormal(obj)
            HalfWidth = obj.ZExtreme * sqrt(obj.NPQ);
            minX = max(0,floor(obj.NP-HalfWidth));
            maxX = min(obj.N,floor(obj.NP+HalfWidth)+1);
            Xs = (minX:maxX);
            CDFs = obj.Standard_Normal.CDF( (Xs+0.5-obj.NP)/sqrt(obj.NPQ) );
            PDFs = diff([0 CDFs]);
        end
        
        function [Xs, PDFs, CDFs] = MakeTablePoisson(obj)
            switch obj.Approx
                case 2
                    obj.PoissonBasis.ResetParms(obj.N*obj.P);
                    Xs   = obj.PoissonBasis.DiscreteX;
                    PDFs = obj.PoissonBasis.DiscretePDF;
                    CDFs = obj.PoissonBasis.DiscreteCDF;
                case 3
                    obj.PoissonBasis.ResetParms(obj.N*obj.Q);
                    Xs   = flip( obj.N - obj.PoissonBasis.DiscreteX ) ;
                    PDFs = flip( obj.PoissonBasis.DiscretePDF );
                    CDFs = 1 - flip( obj.PoissonBasis.DiscreteCDF ) + PDFs;
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
            thisval = unirands<obj.P;
            if obj.N>1
                thissize = size(thisval);
                thisval = squeeze(sum(thisval,numel(thissize)));
            end
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


