classdef Hypergeometric < dDiscrete
    % Hypergeometric(N,K,n): see https://en.wikipedia.org/wiki/Hypergeometric_distribution
    % Consider an urn with N balls of which K are white (others are black).
    % We draw n balls from the urn _without replacement_.
    % For each k=0...min(n,K), what is the probability of drawing exactly k white balls.
    
    properties(SetAccess = protected)
        N, K, n
        MinPDF  % Criterion for deciding when the PDF is small enough to cut off the distribution
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            RealN = NumTrans.GT2Real(1,Parms(1));
            RealK = NumTrans.GT2Real(1,Parms(2));
            Realn = NumTrans.GT2Real(1,Parms(3));
            Reals = [RealN, RealK, Realn];
        end
        
        function Parms = RealsToParms(Reals,~)
            ParmN = NumTrans.GT2Real(1,Reals(1));
            ParmK = NumTrans.GT2Real(1,Reals(2));
            Parmn = NumTrans.GT2Real(1,Reals(3));
            Parms = [ParmN, ParmK, Parmn];
        end
        
    end
    
    methods
        
        function obj=Hypergeometric(varargin)
            obj=obj@dDiscrete('Hypergeometric');
            obj.DistType = 'd';
            obj.ParmTypes = 'iii';
            obj.DefaultParmCodes = 'ffi';
            obj.NDistParms = 3;
            obj.MinPDF = 1e-10;  % was e-20
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Hypergeometric:Constructor', ...
                        'Hypergeometric constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsD(obj);
            obj.K = VerifyIntegerGE(obj,1,newparmvalues(2));
            obj.N = VerifyIntegerGE(obj,obj.K,newparmvalues(1));
            obj.n = VerifyIntegerInRange(obj,0,obj.N,newparmvalues(3));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newN  = ifelse(ParmCodes(1)=='f', obj.N, obj.N+1);
            newK  = ifelse(ParmCodes(2)=='f', obj.K, obj.K+1);
            newn  = ifelse(ParmCodes(3)=='f', obj.n, obj.n+1);
            obj.ResetParms([newN, newK, newn]);
        end
        
        function []=ReInit(obj)
            MakeTables(obj);
            TrimTables(obj,obj.MinPDF,1-eps(1));
            SetBinEdges(obj);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function []=MakeTables(obj)
            obj.LowerBound = max(0,obj.n+obj.K-obj.N);
            obj.UpperBound = min(obj.n,obj.K);
            obj.NValues = round(obj.UpperBound - obj.LowerBound) + 1;
            obj.DiscreteX = (obj.LowerBound:obj.UpperBound)';  % must be a column vector
            obj.DiscretePDF = zeros(obj.NValues,1);
            for ii=1:obj.NValues
                i = obj.DiscreteX(ii);
                obj.DiscretePDF(ii) = nchoosek(obj.K,i) * nchoosek(obj.N-obj.K,obj.n-i);
            end
            obj.DiscretePDF = obj.DiscretePDF / nchoosek(obj.N,obj.n);
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.DiscreteCDF(end) = 1;
        end

        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.n*obj.K/obj.N;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.n*obj.K/obj.N^2*(obj.N-obj.K)*(obj.N-obj.n)/(obj.N-1);
        end
  
%       function thisval=Skewness(obj)
%           if ~obj.Initialized
%               error(UninitializedError(obj));
%           end
%           thisval = SEE WIKIPEDIA
%       end
  
%       function thisval=Kurtosis(obj)
%           if ~obj.Initialized
%               error(UninitializedError(obj));
%           end
%           thisval = 9 + SEE WIKIPEDIA
%       end
  
    end  % methods
    
end  % class Hypergeometric



