classdef LinearTrans < dTransMono
    % LinearTrans(BasisRV,Multiplier,Addend) produces the linear transformation BasisRV*Multiplier+Addend.
    
    properties(SetAccess = protected)
        Multiplier, Addend
    end
    
    methods (Static)
        
        function TransReals = TransParmsToReals(Parms,~)
            TransReals = Parms(end-1:end);
        end
        
        function TransParms = TransRealsToParms(Reals,~)
            TransParms = Reals(end-1:end);
        end
        
    end
    
    methods
        
        function obj=LinearTrans(BasisDist,Multiplier,Addend)
            obj=obj@dTransMono('LinearTrans',BasisDist);
            obj.Multiplier = Multiplier;
            obj.Addend = Addend;
            obj.TransReverses = Multiplier < 0;
            obj.PDFScaleFactorKnown = true;
            obj.AddParms(2,'rr');
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Multiplier = newparmvalues(end-1);
            obj.Addend = newparmvalues(end);
            ResetParms@dTransMono(obj,newparmvalues(1:end-2));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewMultiplier = ifelse(ParmCodes(end-1)=='f', obj.Multiplier,0.95*obj.Multiplier);
            NewAddend = ifelse(ParmCodes(end)=='f', obj.Addend,obj.Addend+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewMultiplier NewAddend]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.Multiplier obj.Addend];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX * obj.Multiplier + obj.Addend;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = (TransX - obj.Addend) / obj.Multiplier;
        end
        
        function thisval = PDFScaleFactor(obj,~)
            thisval = abs(1/obj.Multiplier);
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Mean(obj.BasisRV) * obj.Multiplier + obj.Addend;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Variance(obj.BasisRV) * obj.Multiplier^2;
        end
        
    end  % methods
    
end  % class LinearTrans


