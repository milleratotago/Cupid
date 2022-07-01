classdef TruncatedP < TruncParent
    % TruncatedP(BasisRV,LowerCutoffP,UpperCutoffP) produces the truncated version of the BasisRV,
    % truncating between the two indicated CDF cutoffs.
    %
    % You must specify fixed cutoffs in the constructor if you want them, like this:
    %   truncdist = TruncatedP(dist,.01,.98,'FixedCutoffLow','FixedCutoffHi');
    
    methods
        
        function obj=TruncatedP(BasisDist,LowerP,UpperP,varargin)
            obj=obj@TruncParent('TruncatedP',BasisDist,varargin{:});
            obj.NewCutoffs(obj.BasisRV.InverseCDF(LowerP),obj.BasisRV.InverseCDF(UpperP));
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
            if obj.FixedCutoffHi
                HiP = obj.UpperCutoffP;
            else
                HiP = newparmvalues(end);
                newparmvalues = newparmvalues(1:end-1);  % remove it so now LowP is on the end
            end
             if obj.FixedCutoffLow
                LowP = obj.LowerCutoffP;
            else
                LowP = newparmvalues(end);
            end
            obj.NewCutoffs(obj.BasisRV.InverseCDF(LowP),obj.BasisRV.InverseCDF(HiP));
            obj.Initialized = true;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            tempParms = obj.BasisRV.ParmValues;
            if ~obj.FixedCutoffLow
                NewLowerP = ifelse(ParmCodes(end-1)=='f', obj.LowerCutoffP,0.9*obj.LowerCutoffP);
                tempParms = [tempParms NewLowerP];
            end
            if ~obj.FixedCutoffHi
                NewUpperP = ifelse(ParmCodes(end)=='f', obj.UpperCutoffP,obj.UpperCutoffP+0.1*(1-obj.UpperCutoffP));
                tempParms = [tempParms NewUpperP];
            end
            obj.ResetParms(tempParms);
        end
        
        function parmvals = TransParmValues(obj)
            if ~obj.FixedCutoffLow
                parmvals = obj.LowerCutoffP;
            else
                parmvals = [];
            end
            if ~obj.FixedCutoffHi
                parmvals = [parmvals obj.UpperCutoffP];
            end
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            % TransReals = NumTrans.Bounded2Real(0,1,Parms(end-1:end));
            if ~obj.FixedCutoffLow && ~obj.FixedCutoffHi
                TransReals(2) = NumTrans.Bounded2Real(0,1,Parms(end));
                TransReals(1) = NumTrans.Bounded2Real(0,Parms(end),Parms(end-1));
            elseif obj.FixedCutoffLow && obj.FixedCutoffHi
                TransReals = [];
            elseif obj.FixedCutoffLow
                TransReals = NumTrans.Bounded2Real(0,1,Parms(end));
            else
                TransReals = NumTrans.Bounded2Real(0,obj.UpperCutoffP,Parms(end));
            end
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            % TransParms = NumTrans.Real2Bounded(0,1,Reals(end-1:end));
            if ~obj.FixedCutoffLow && ~obj.FixedCutoffHi
                TransParms(2) = NumTrans.Real2Bounded(0,1,Reals(end));
                TransParms(1) = NumTrans.Real2Bounded(0,TransParms(2),Reals(end-1));
            elseif obj.FixedCutoffLow && obj.FixedCutoffHi
                TransParms = [];
            elseif obj.FixedCutoffLow
                TransParms = NumTrans.Real2Bounded(0,1,Reals(end));
            else
                TransParms = NumTrans.Real2Bounded(0,obj.UpperCutoffP,Reals(end));
            end
        end

        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' obj.BasisRV.StringName];
            obj.StringName = [obj.StringName ',' num2str(obj.LowerCutoffP) ',' num2str(obj.UpperCutoffP)];
            obj.StringName = [obj.StringName ')'];
        end
        
    end  % methods
    
end  % class TruncatedP



