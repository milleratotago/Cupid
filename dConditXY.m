classdef dConditXY < dContinuous
    % An abstract class for examining the distribution of random variable X conditional
    % on its being either less than or greater than an independent random variable Y

    properties(SetAccess = protected)
        BasisRV1, BasisRV2
    end

    methods(Abstract)
        [] = ReInit(obj)
        thispdf = PDF(obj,X)
    end

    methods

        function obj=dConditXY(FamilyName, BasisRV1, BasisRV2)   % Constructor
            obj=obj@dContinuous(FamilyName);  % Inherited constructor
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
