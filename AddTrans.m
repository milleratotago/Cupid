classdef AddTrans < dTransMono
    % AddTrans(BasisRV,Addend): Shifts the BasisRV by the specified additive constant.
    
    properties(SetAccess = protected)
        Addend
    end
    
    methods
        
        function obj=AddTrans(BasisDist,Addend)
            obj=obj@dTransMono('AddTrans',BasisDist);
            obj.Addend = Addend;
            obj.AddParms(1,'r');
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Addend = newparmvalues(end);
            ResetParms@dTransMono(obj,newparmvalues(1:end-1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewAddend = ifelse(ParmCodes(end)=='f', obj.Addend,obj.Addend+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewAddend]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.Addend;
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX + obj.Addend;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX - obj.Addend;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
        function thisval = PDFScaleFactor(~,~)
            thisval = 1;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Mean(obj.BasisRV) + obj.Addend;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Variance(obj.BasisRV);
        end

        function thisval=Skewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Skewness(obj.BasisRV);
        end

        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Kurtosis(obj.BasisRV);
        end

        function thisval=CenMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = CenMoment(obj.BasisRV,I);
        end

    end  % methods
    
end  % class AddTrans


