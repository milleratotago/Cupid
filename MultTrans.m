classdef MultTrans < dTransMono
    % MultTrans(BasisRV,Multiplier): Rescales the BasisRV by a constant multiplier.
    
    properties(SetAccess = protected)
        Multiplier
    end
    
    methods
        
        function obj=MultTrans(BasisDist,Multiplier)
            obj=obj@dTransMono('MultTrans',BasisDist);
            obj.Multiplier = Multiplier;
            obj.AddParms(1,'r');
            obj.TransReverses = Multiplier < 0;
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Multiplier = newparmvalues(end);
            ResetParms@dTransMono(obj,newparmvalues(1:end-1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewMult = ifelse(ParmCodes(end)=='f', obj.Multiplier,obj.Multiplier*1.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewMult]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.Multiplier;
        end
        
        % Making the next 2 static causes problems:
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX * obj.Multiplier;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX / obj.Multiplier;
        end
        
        function thisval = PDFScaleFactor(obj,~)
            thisval = 1/obj.Multiplier;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Mean(obj.BasisRV) * obj.Multiplier;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Variance(obj.BasisRV) * obj.Multiplier^2;
        end
        
    end  % methods
    
end  % class MultTrans


