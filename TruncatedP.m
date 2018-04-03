classdef TruncatedP < TruncatedX
    % TruncatedP(BasisRV,LowerCutoffP,UpperCutoffP) produces the truncated version of the BasisRV,
    % truncating between the two indicated CDF cutoffs.
    
    methods
        
        function obj=TruncatedP(varargin)
            obj=obj@TruncatedX;
            obj.ThisFamilyName = 'TruncatedP';
            switch nargin
                case 0
                case 3
                    BuildMyBasis(obj,varargin{1});
                    ResetParms(obj,[ParmValues(obj.BasisRV) varargin{end-1} varargin{end}]);
                otherwise
                    ME = MException('TruncatedP:Constructor', ...
                        'TruncatedP constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.ThisFamilyName '(' obj.BasisRV.StringName ',' num2str(obj.LowerCutoffP)  ',' num2str(obj.UpperCutoffP) ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues(1:end-2));
            obj.LowerCutoffP = newparmvalues(end-1);
            obj.UpperCutoffP = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLowerP = ifelse(ParmCodes(end-1)=='f', obj.LowerCutoffP,0.9*obj.LowerCutoffP);
            NewUpperP = ifelse(ParmCodes(end)=='f', obj.UpperCutoffP,obj.UpperCutoffP+0.1*(1-obj.UpperCutoffP));
            obj.ResetParms([obj.BasisRV.ParmValues NewLowerP NewUpperP]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.UnconditionalP = obj.UpperCutoffP - obj.LowerCutoffP;
            assert(obj.UnconditionalP>0,'TruncatedP distribution must include probability > 0.');
            assert(obj.BasisRV.Initialized,'TruncatedP BasisRV must already be initialized.');
            obj.LowerCutoffX = obj.BasisRV.InverseCDF(obj.LowerCutoffP);
            obj.UpperCutoffX = obj.BasisRV.InverseCDF(obj.UpperCutoffP);
            ReInit@TruncatedX(obj);
            % if (obj.NameBuilding)
            %     BuildMyName(obj);
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.BasisRV.ParmValues obj.LowerCutoffP obj.UpperCutoffP];
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



