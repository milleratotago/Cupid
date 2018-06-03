classdef TruncatedP < TruncParent
    % TruncatedP(BasisRV,LowerCutoffP,UpperCutoffP) produces the truncated version of the BasisRV,
    % truncating between the two indicated CDF cutoffs.
    
    methods
        
        function obj=TruncatedP(BasisDist,LowerP,UpperP)
            obj=obj@TruncParent('TruncatedP',BasisDist);
            obj.NewCutoffs(obj.BasisRV.InverseCDF(LowerP),obj.BasisRV.InverseCDF(UpperP));
            obj.AddParms(2,'rr');
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NewCutoffs(obj.BasisRV.InverseCDF(newparmvalues(end-1)),obj.BasisRV.InverseCDF(newparmvalues(end)));
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLowerP = ifelse(ParmCodes(end-1)=='f', obj.LowerCutoffP,0.9*obj.LowerCutoffP);
            NewUpperP = ifelse(ParmCodes(end)=='f', obj.UpperCutoffP,obj.UpperCutoffP+0.1*(1-obj.UpperCutoffP));
            obj.ResetParms([obj.BasisRV.ParmValues NewLowerP NewUpperP]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.LowerCutoffP obj.UpperCutoffP];
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            % TransReals = NumTrans.Bounded2Real(0,1,Parms(end-1:end));
            TransReals(2) = NumTrans.Bounded2Real(0,1,Parms(end));
            TransReals(1) = NumTrans.Bounded2Real(0,Parms(end),Parms(end-1));
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            % TransParms = NumTrans.Real2Bounded(0,1,Reals(end-1:end));
            TransParms(2) = NumTrans.Real2Bounded(0,1,Reals(end));
            TransParms(1) = NumTrans.Real2Bounded(0,TransParms(2),Reals(end-1));
        end
        
    end  % methods
    
end  % class TruncatedP



