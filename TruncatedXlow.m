classdef TruncatedXlow < TruncParent
    % TruncatedXlow(BasisRV,LowerCutoffX) produces the truncated version of the BasisRV,
    % truncating at the low end with the indicated X cutoff.
    
    methods
        
        function obj=TruncatedXlow(BasisDist,LowerX,varargin)
            obj=obj@TruncParent('TruncatedXlow',BasisDist,varargin{:});
            obj.NewCutoffs(LowerX,obj.BasisRV.UpperBound);
            if ~obj.FixedCutoffLow
                obj.AddParms(1,'r');
            end
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NewCutoffs(newparmvalues(end),obj.BasisRV.UpperBound);
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLowerX = ifelse(ParmCodes(end)=='f', obj.LowerCutoffX,obj.LowerCutoffX-0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewLowerX]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.LowerCutoffX;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            if ~obj.FixedCutoffLow
                TransReals = Parms(end);
            else
                TransReals = [];
            end
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            if ~obj.FixedCutoffLow
                TransParms = Reals(end);
            else
                TransParms = [];
            end
        end
        
    end  % methods
    
end  % class TruncatedXlow



