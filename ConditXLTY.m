classdef ConditXLTY < dContinuous
    % Distribution of random variable X conditional on its being less than an independent random variable Y

    properties(SetAccess = protected)
        BasisRV1, BasisRV2
        PrXLTY
    end

    methods

        function obj=ConditXLTY(BasisRV1, BasisRV2)   % Constructor
            obj=obj@dContinuous('ConditXLTY');  % Inherited constructor
            obj.BasisRV1 = BasisRV1;
            obj.BasisRV2 = BasisRV2;
            assert( (BasisRV1.DistType=='c') && (BasisRV2.DistType=='c'), ...
              [obj.FamilyName ' can only handle continuous distributions (so far).']);
            obj.DistType = BasisRV1.DistType;
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
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.PrXLTY = 1 - PrXGTY(obj.BasisRV1,obj.BasisRV2);
            obj.LowerBound = obj.BasisRV1.LowerBound;
            obj.UpperBound = min(obj.BasisRV1.UpperBound,obj.BasisRV2.UpperBound);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.BasisRV1.PDF(X(InBounds)) .* (1 - obj.BasisRV2.CDF(X(InBounds))) / obj.PrXLTY;
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
        
    end

end % classdef
