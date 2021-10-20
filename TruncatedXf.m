classdef TruncatedXf < TruncParent  % NEWJEFF: NEW, NO UNITTESTS OR DOCS This version has fixed cutoffs--not adjustable parameters.
    % TruncatedXf(BasisRV,LowerCutoffX,UpperCutoffX) produces the truncated version of the BasisRV,
    % truncating between the two indicated X cutoffs.
    
    methods
        
        function obj=TruncatedXf(BasisDist,LowerX,UpperX)
            obj=obj@TruncParent('TruncatedXf',BasisDist);
            obj.NewCutoffs(LowerX,UpperX);
            % obj.AddParms(2,'rr');
            obj.ReInit;
        end

        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            % obj.NewCutoffs(newparmvalues(end-1),newparmvalues(end));
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes);
            obj.BasisRV.PerturbParms(ParmCodes);
            % NewLowerX = ifelse(ParmCodes(end-1)=='f', obj.LowerCutoffX,obj.LowerCutoffX-0.1);
            % NewUpperX = ifelse(ParmCodes(end)=='f', obj.UpperCutoffX,obj.UpperCutoffX+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues obj.LowerCutoffX obj.UpperCutoffX]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [];
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = [];
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = [];
        end
        
    end  % methods
    
end  % class TruncatedXf



