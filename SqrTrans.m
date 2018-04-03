classdef SqrTrans < PowerTrans
    % SqrTrans(BasisRV): Sqr transformation of a BasisRV which must be NONNEGATIVE.
    
    methods
        
        function obj=SqrTrans(varargin)
            assert(nargin==1,'SqrTrans constructor needs one argument (a distribution).');
            obj=obj@PowerTrans(varargin{:},2);
            obj.ThisFamilyName = 'SqrTrans';
            obj.NTransParms = 0;
            obj.TransParmCodes = '';
            obj.DistType = obj.BasisRV.DistType;
            obj.NDistParms = obj.BasisRV.NDistParms;
            obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
            obj.BuildMyName;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@PowerTrans(obj,[newparmvalues obj.Power]);
            ReInit(obj);
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.ThisFamilyName '(' obj.BasisRV.StringName ')'];
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            obj.ResetParms(obj.BasisRV.ParmValues);
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
    
end  % class SqrTrans



