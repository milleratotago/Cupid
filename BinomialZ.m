classdef BinomialZ < dDiscrete  % Do not descend from Binomial.
    % BinomialZ(N,P,PerfectAdjust): Transformation of a binomial(N,P) where each binomial K
    % is converted to the Z score with cumulative probability k/N.
    % PerfectAdjust indicates what to do if k=0 or k=N: replace 0 with PerfectAdjust and N with N-PerfectAdjust.
    
    properties
        N
        P
        PerfectAdjust
        myBinomial
    end
    
    methods
        
        function obj=BinomialZ(varargin)
            obj=obj@dDiscrete('BinomialZ');
            obj.N = varargin{1};
            obj.P = varargin{2};
            obj.PerfectAdjust = varargin{3};
            obj.myBinomial = Binomial(obj.N,obj.P);
            obj.ParmTypes = [obj.myBinomial.ParmTypes 'r'];
            obj.DefaultParmCodes = [obj.myBinomial.DefaultParmCodes 'f'];
            obj.NDistParms = obj.myBinomial.NDistParms + 1;
            BuildMyName(obj);
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsD(obj);
            obj.N = round(newparmvalues(1));
            obj.P = newparmvalues(2);
            obj.PerfectAdjust = newparmvalues(3);
            obj.myBinomial.ResetParms([obj.N,obj.P]);
            ReInit(obj);
        end
        
        function []=PerturbParms(obj,ParmCodes)
            obj.myBinomial.PerturbParms(ParmCodes);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            obj.myBinomial.ReInit;
            obj.MakeTables;
            if obj.NameBuilding
                BuildMyName(obj);
            end
        end
        
        function []=MakeTables(obj)
            obj.NValues = obj.myBinomial.NValues;
            obj.DiscretePDF = obj.myBinomial.DiscretePDF;
            obj.DiscreteCDF = obj.myBinomial.DiscreteCDF;
            obj.DiscreteCDFmax = obj.myBinomial.DiscreteCDFmax;
            obj.DiscreteX = obj.myBinomial.DiscreteX;
            if obj.DiscreteX(1) == 0
                obj.DiscreteX(1) = obj.PerfectAdjust;
            end
            if obj.DiscreteX(end) == obj.myBinomial.N
                obj.DiscreteX(end) = obj.myBinomial.N - obj.PerfectAdjust;
            end
            obj.DiscreteX = norminv(obj.DiscreteX / obj.N);
%           obj.SetBinEdges;  % The next 3 lines work better.
            obj.DiscreteXmin = obj.DiscreteX - 100*eps(obj.DiscreteX);  % NEWJEFF: Use XGrain here
            obj.DiscreteXmax = obj.DiscreteX + 100*eps(obj.DiscreteX);
            obj.StoredTablesInitialized = true;
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
            obj.Initialized = true;
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.myBinomial.ParmsToReals(Parms(1:2)) NumTrans.Bounded2Real(0,1,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.myBinomial.RealsToParms(Reals(1:2)) NumTrans.Real2Bounded(0,1,Reals(3))];
        end
        
    end  % methods
    
end  % class BinomialZ


