classdef SqrtTrans < PowerTrans
    % SqrtTrans(BasisRV): Sqrt transformation of a BasisRV which must be NONNEGATIVE.
    
    methods
        
        function obj=SqrtTrans(varargin)
            assert(nargin==1,'SqrtTrans constructor needs one argument (a distribution).');
            obj=obj@PowerTrans(varargin{:},0.5);
            obj.ThisFamilyName = 'SqrtTrans';
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
    
end  % class SqrtTrans



