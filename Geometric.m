classdef Geometric < dDiscrete     % NWJEFF: Not vectorized
    % Geometric(P): Number of trials to first success, including the success (i.e., X=1,2,3,...)
    
    properties(SetAccess = protected)
        P, Q
        MinPDF  % Criterion for deciding when the PDF is small enough to cut off the distribution
    end
    
    methods
        
        function obj=Geometric(varargin)
            obj=obj@dDiscrete('Geometric');
            obj.DistType = 'd';
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            obj.MinPDF = 1e-10;  % was e-20
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('Geometric:Constructor', ...
                        'Geometric constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsD(obj);
            obj.P = newparmvalues(1);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newP  = ifelse(ParmCodes(1)=='f', obj.P, obj.P*0.95);
            obj.ResetParms(newP);
        end
        
        function []=ReInit(obj)
            assert(obj.P>0&&obj.P<1,'Geometric P must be 0-1.');
            MakeTables(obj);
            TrimTables(obj,obj.MinPDF,1-eps(1));
            SetBinEdges(obj);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function []=MakeTables(obj)
            obj.Q = 1 - obj.P;
            obj.LowerBound = 1;
            obj.UpperBound = ceil(1+log(obj.MinPDF/obj.P)/log(1-obj.P));
            obj.NValues = round(obj.UpperBound - obj.LowerBound) + 1;
            obj.DiscreteX = obj.LowerBound:obj.UpperBound;
            AllQs = obj.Q*ones(1,obj.NValues);
            obj.DiscretePDF = obj.P*cumprod(AllQs)/obj.Q;
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.DiscreteCDF(end) = 1;
        end

        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.Bounded2Real(0,1,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2Bounded(0,1,Reals(1));
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1/obj.P;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (1-obj.P)/obj.P^2;
        end
  
        function thisval=Skewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (2-obj.P)/sqrt(1-obj.P);
        end
  
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 9 + obj.P^2/(1-obj.P);
        end
  
        function s=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            meanObs = mean(Observations);
            ResetParms(obj,1/meanObs);
            BuildMyName(obj);
            s=obj.StringName;
        end
        
        function thisval=Random(obj,varargin)
            thisval = geornd(obj.P,varargin{:}) + 1;  % geornd is number of failures before first success
        end
        
    end  % methods
    
end  % class Geometric



