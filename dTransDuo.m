classdef dTransDuo < dEither   % NWJEFF: Need a demo of using this generically for researcher's own novel function.
    % dTransDuo(BasisRV1,BasisRV2) is an abstract class to create a random variable that is
    % a function of two independent basis random variables, BasisRV1+BasisRV2.
    
    properties(SetAccess = protected)
        BasisRV1, BasisRV2
        CDFRelTol, CDFAbsTol
    end

    methods (Abstract)

      FNXY = FofDuo(X,Y)  % e.g., X+Y for Convolution
      X = ReverseFofDuo(FNXY,Y)  % e.g., FNXY-Y for Convolution

    end

    
    methods
        
        function obj=dTransDuo(FamilyName,Basis1,Basis2)
            obj=obj@dEither(FamilyName);
            obj.BasisRV1 = Basis1;
            obj.BasisRV2 = Basis2;
            obj.CDFRelTol = 1e-3;
            obj.CDFAbsTol = 1e-4;
            assert(Basis1.DistType==Basis2.DistType, ...
              [obj.FamilyName ' can only handle Basis distributions of the same type (so far).']);
            obj.DistType = Basis1.DistType;
            obj.NDistParms = obj.BasisRV1.NDistParms + obj.BasisRV2.NDistParms;
            obj.DefaultParmCodes = [obj.BasisRV1.DefaultParmCodes obj.BasisRV2.DefaultParmCodes];
            obj.ReInit;
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV1.PerturbParms(ParmCodes(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.PerturbParms(ParmCodes(obj.BasisRV1.NDistParms+1:end));
            obj.ResetParms([obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
        end
        
        function BuildMyName(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            obj.StringName = [obj.FamilyName '(' obj.BasisRV1.StringName ',' obj.BasisRV2.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj);
            obj.BasisRV1.ResetParms(newparmvalues(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.ResetParms(newparmvalues(obj.BasisRV1.NDistParms+1:end));
            ReInit(obj);
        end
        
        function MakeTables(obj)
            [Pairs, Probs] = PairsNProbs(obj.BasisRV1.DiscreteX,obj.BasisRV2.DiscreteX,obj.BasisRV1.DiscretePDF,obj.BasisRV2.DiscretePDF);
            Results = FofDuo(obj,Pairs(1,:),Pairs(2,:));
            [obj.DiscreteX, obj.DiscretePDF] = CollapseVals(Results,Probs);
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.DiscreteCDF(end) = 1;
            obj.NValues = numel(obj.DiscreteX);
            obj.SetBinEdges;
            obj.LowerBound = obj.DiscreteXmin(1);
            obj.UpperBound = obj.DiscreteXmax(end);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            switch obj.DistType
                case 'd'
                    obj.MakeTables;
                case 'c'
                    obj.SetBoundsContin;
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thiscdf=CDF(obj,X)
            if obj.DistType=='d'
                thiscdf = CDF@dDiscrete(obj,X);
                return;
            end
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    fun = @(x) obj.BasisRV1.CDF(obj.ReverseFofDuo(X(iel),x)).*obj.BasisRV2.PDF(x);
                    thiscdf(iel)=integral(fun,obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,'RelTol', obj.CDFRelTol, 'AbsTol', obj.CDFAbsTol);
                end
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.BasisRV1.ParmsToReals(Parms(1:obj.BasisRV1.NDistParms)) obj.BasisRV2.ParmsToReals(Parms(obj.BasisRV1.NDistParms+1:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.BasisRV1.RealsToParms(Reals(1:obj.BasisRV1.NDistParms)) obj.BasisRV2.RealsToParms(Reals(obj.BasisRV1.NDistParms+1:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues];
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.FofDuo(obj.BasisRV1.Random(varargin{:}),obj.BasisRV2.Random(varargin{:}));
        end
        
    end  % methods
    
end  % class dTransDuo

