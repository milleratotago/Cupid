classdef TruncatedXlow < TruncParent
    % TruncatedXlow(BasisRV,LowerCutoffX) produces the truncated version of the BasisRV,
    % truncating at the low end with the indicated X cutoff.
    
    methods
        
        function obj=TruncatedXlow(BasisDist,LowerX)
            obj=obj@TruncParent('TruncatedXlow',BasisDist);
            obj.NewCutoffs(LowerX,obj.BasisRV.UpperBound);
            obj.AddParms(1,'r');
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NewCutoffs(newparmvalues(end),obj.BasisRV.UpperBound);
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes);
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLowerX = ifelse(ParmCodes(end)=='f', obj.LowerCutoffX,obj.LowerCutoffX-0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewLowerX]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.LowerCutoffX;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
    end  % methods
    
end  % class TruncatedXlow



