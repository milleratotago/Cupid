classdef BinomialMixed < dDiscrete  % NEWJEFF: Undocumented distribution.
    % BinomialMixed(P): Distribution of the number of successes in N trials,
    %   where the probability of success varies across trials.
    % P is a vector of length (1,N) showing the Pr_of_success for each of the N trials.
    % N cannot be varied.
    
    % Useful reference: https://stats.stackexchange.com/questions/9510/probability-distribution-for-different-probabilities
    
    properties(SetAccess = protected)
        N, P
    end
    
    methods
        
        function obj=BinomialMixed(varargin)
            obj=obj@dDiscrete('BinomialMixed');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'f';  % NWJEFF: Probability adjustment not supported
            obj.NDistParms = 1;
            obj.UseStoredTables = true;
            obj.Grain = 10;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('BinomialMixed:Constructor', ...
                        'BinomialMixed constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.P = newparmvalues(:);
            ReInit(obj);
        end
        
        function []=PerturbParms(obj,ParmCodes)
            obj.P = ifelse(ParmCodes(1)=='f',obj.P,0.95*obj.P);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            B = obj.P(obj.P<0 | obj.P>1);
            assert(numel(B)==0,'BinomialMixed P must be in range (0-1).');
            obj.N = numel(obj.P);
            obj.LowerBound = 0;
            obj.UpperBound =  obj.N;
            obj.Initialized = true;
            obj.LowerBound = 0;
            obj.UpperBound = obj.N;
            obj.MakeTables;
            i=obj.LowerBound;
            while obj.CDF(i)<=eps(0)
                i=i+1;
            end
            obj.LowerBound = max(obj.LowerBound,i-1);
            i=obj.UpperBound;
            while obj.CDF(i)>=1-eps
                i=i-1;
            end
            obj.UpperBound = min(obj.UpperBound,i+1);
            obj.NValues = obj.UpperBound - obj.LowerBound + 1;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.Bounded2Real(0,1,Parms(1))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2Bounded(0,1,Reals(1))];
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
        
        function thisval=NextValue(obs,X)
            thisval = X+1;
        end
        
        function []=MakeTables(obj)
            obj.NValues = obj.N + 1;  % first element is zero
            obj.StoredX = 0:obj.N;
            obj.StoredPDF = obj.StoredX;
            obj.StoredCDF = obj.StoredX;
            
            syms tempt;
            pgfi = 1-obj.P*(1-tempt);
            pgfS = collect( prod(pgfi) );  % collecting terms is essential!
            
            RunningTotal = 0;
            for i=0:obj.N
                if i==0
                    dpgfS = pgfS;
                    ifac = 1;
                else
                    dpgfS = diff(dpgfS);
                    % dpgfS = collect(dpgfS,tempt); % Seems unnecessary
                    ifac = ifac*i;
                end
                obj.StoredPDF(i+1) = double( 1/ifac * limit(dpgfS,tempt,0) );
                RunningTotal = RunningTotal + obj.StoredPDF(i+1);
                obj.StoredCDF(i+1) = RunningTotal;
            end
            obj.StoredTablesInitialized = true;
        end
        
        function thisval=Random(obj,varargin)
            if numel(varargin)==0
                varargin{1} = 1;
            end
            unirands = rand(obj.N,varargin{:});
            thisvar01 = unirands<obj.P;
            thisval = squeeze(sum(thisvar01));  % NWJEFF: Transposes rows and columns
        end
        
    end  % methods
    
end  % class BinomialMixed


