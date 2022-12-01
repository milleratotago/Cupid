classdef TruncatedX < TruncParent
    % TruncatedX(BasisRV,LowerCutoffX,UpperCutoffX) produces the truncated version of the BasisRV,
    % truncating between the two indicated X cutoffs.
    
    methods
        
        function obj=TruncatedX(BasisDist,LowerX,UpperX,varargin)
            obj=obj@TruncParent('TruncatedX',BasisDist,varargin{:});
            obj.LowerCutoffX = LowerX;
            obj.UpperCutoffX = UpperX;
            obj.NewCutoffXs(LowerX,UpperX);
            if ~obj.FixedCutoffLow
                obj.AddParms(1,'r');
            end
            if ~obj.FixedCutoffHi
                obj.AddParms(1,'r');
            end
            obj.ReInit;
        end

        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NewCutoffXs(newparmvalues(end-1),newparmvalues(end));
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes);
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLowerX = ifelse(ParmCodes(end-1)=='f', obj.LowerCutoffX,obj.LowerCutoffX-0.1);
            NewUpperX = ifelse(ParmCodes(end)=='f', obj.UpperCutoffX,obj.UpperCutoffX+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewLowerX NewUpperX]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.LowerCutoffX obj.UpperCutoffX];
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            if ~obj.FixedCutoffLow && ~obj.FixedCutoffHi
                TransReals = Parms(end-1:end);
            elseif obj.FixedCutoffLow && obj.FixedCutoffHi
                TransReals = [];
            else
                TransReals = Parms(end);  % Either low or hi
            end
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            if ~obj.FixedCutoffLow && ~obj.FixedCutoffHi
               TransParms = Reals(end-1:end);
            elseif obj.FixedCutoffLow && obj.FixedCutoffHi
                TransParms = [];
            else
                TransParms = Reals(end);  % Either low or hi
            end
        end
        
    end  % methods
    
end  % class TruncatedX



