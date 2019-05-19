classdef BinomialMixed < dDiscrete
    % BinomialMixed(P): Distribution of the number of successes in N trials,
    %   where the probability of success varies across trials.
    % P is a vector of length (1,N) showing the Pr_of_success for each of the N trials.
    % N cannot be varied.
    
    % Useful reference: https://stats.stackexchange.com/questions/9510/probability-distribution-for-different-probabilities
    
    properties(SetAccess = protected)
        N, P
    end

    properties(SetAccess = public)
        NumDP  % N of decimal places when converting probabilities to string in distribution name.
    end

    methods
        
        function obj=BinomialMixed(varargin)
            obj=obj@dDiscrete('BinomialMixed');
            obj.NumDP = 3;
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
        
        function BuildMyName(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            s = mat2str(obj.P,obj.NumDP);
            s = strrep(s,' ',',');
            obj.StringName = ['BinomialMixed(' s ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            if isequal(newparmvalues,obj.P)
                return;  % Do nothing if nothing has changed.
            end
            ClearBeforeResetParmsD(obj);
            obj.N = numel(newparmvalues);
            obj.P = zeros(1,obj.N);  % Make sure it is a row vector.
            obj.P(:) = newparmvalues(:);
            obj.NDistParms = obj.N;  % The vector of p's is not adjustable, but we need some parameters.
            obj.ParmTypes = repmat('r',1,obj.N);
            obj.DefaultParmCodes = repmat('f',1,obj.N);  % NWJEFF: Probability adjustment not supported, so this distribution has no parameters
            ReInit(obj);
        end
        
        function parmvals = ParmValues(obj)
            parmvals = obj.P;  % Not adjustable so don't return it as a parameter.
        end
        
        function []=PerturbParms(obj,ParmCodes)
            for i=1:obj.N
                obj.P(i) = ifelse(ParmCodes(i)=='f',obj.P(i),0.95*obj.P(i));
            end
            if numel(strfind(ParmCodes,'r'))>0
                ReInit(obj);  % Only if something changed.
            end
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            B = obj.P(obj.P<0 | obj.P>1);
            assert(numel(B)==0,'BinomialMixed P must be in range (0-1).');
            % obj.LowerBound = 0;
            % obj.UpperBound =  obj.N;
            obj.MakeTables;
            obj.TrimTables(obj.PDFNearlyZero,1);
            obj.SetBinEdges;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = Parms; % [NumTrans.Bounded2Real(0,1,Parms(1))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = Reals; % [NumTrans.Real2Bounded(0,1,Reals(1))];
        end
        
        function []=MakeTables(obj)
            obj.N = numel(obj.P);
            obj.LowerBound = 0;
            obj.UpperBound = obj.N;
            obj.NValues = obj.N + 1;  % first element is zero
            obj.DiscreteX = 0:obj.N;
            obj.DiscretePDF = obj.DiscreteX;
            obj.DiscreteCDF = obj.DiscreteX;
            
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
                obj.DiscretePDF(i+1) = double( 1/ifac * limit(dpgfS,tempt,0) );
                RunningTotal = RunningTotal + obj.DiscretePDF(i+1);
                obj.DiscreteCDF(i+1) = RunningTotal;
            end
            obj.StoredTablesInitialized = true;
            obj.DiscreteCDF(end) = 1;
        end
        
        function thisval=Random(obj,varargin)
            if numel(varargin)==0
                varargin{1} = 1;
            end
            unirands = rand(obj.N,varargin{:});
            thisvar01 = zeros(size(unirands));
            for i=1:obj.N
                thisvar01(i,:) = unirands(i,:)<=obj.P(i);
            end
%             thisvar01 = unirands<obj.P;
            thisval = squeeze(sum(thisvar01));  % NWJEFF: Transposes rows and columns
        end
        
    end  % methods
    
end  % class BinomialMixed


