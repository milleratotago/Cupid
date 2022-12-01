classdef TruncatedXhi < TruncParent
    % TruncatedXhi(BasisRV,UpperCutoffX) produces the truncated version of the BasisRV,
    % truncating at the high end with the indicated X cutoff.
    
    methods
        
        function obj=TruncatedXhi(BasisDist,UpperX,varargin)
            obj=obj@TruncParent('TruncatedXhi',BasisDist,varargin{:});
            obj.NewCutoffXs(obj.BasisRV.LowerBound,UpperX);
            if ~obj.FixedCutoffHi
                obj.AddParms(1,'r');
            end
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NewCutoffXs(obj.BasisRV.LowerBound,newparmvalues(end));
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
            if ~obj.FixedCutoffHi
                TransReals = Parms(end);
            else
                TransReals = [];
            end
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            if ~obj.FixedCutoffHi
                TransParms = Reals(end);
            else
                TransParms = [];
            end
        end
        
    end  % methods
    
end  % class TruncatedXhi



