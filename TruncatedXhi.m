classdef TruncatedXhi < TruncParent
    % TruncatedXhi(BasisRV,UpperCutoffX) produces the truncated version of the BasisRV,
    % truncating at the high end with the indicated X cutoff.
    
    methods
        
        function obj=TruncatedXhi(BasisDist,UpperX)
            obj=obj@TruncParent('TruncatedXhi',BasisDist);
            obj.NewCutoffs(obj.BasisRV.LowerBound,UpperX);
            obj.AddParms(1,'r');
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NewCutoffs(obj.BasisRV.LowerBound,newparmvalues(end));
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes);
            obj.BasisRV.PerturbParms(ParmCodes);
            NewUpperX = ifelse(ParmCodes(end)=='f', obj.UpperCutoffX,obj.UpperCutoffX+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewUpperX]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.UpperCutoffX;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
    end  % methods
    
end  % class TruncatedXhi



