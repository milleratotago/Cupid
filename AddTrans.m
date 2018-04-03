classdef AddTrans < dTransOf1
    % AddTrans(BasisRV,Addend): Shifts the BasisRV by the specified additive constant.
    
    properties(SetAccess = protected)
        Addend
    end
    
    methods
        
        function obj=AddTrans(varargin)
            obj@dTransOf1('AddTrans');
            obj.NTransParms = 1;
            obj.TransParmCodes = 'r';
            switch nargin
                case 0
                case 2
                    BuildMyBasis(obj,varargin{1});
                    obj.DistType = obj.BasisRV.DistType;
                    obj.NDistParms = obj.BasisRV.NDistParms + 1;
                    obj.DefaultParmCodes = [obj.BasisRV.DefaultParmCodes 'r'];
                    ResetParms(obj,[obj.BasisRV.ParmValues varargin{2}]);
                otherwise
                    ME = MException('AddTrans:Constructor', ...
                        'AddTrans constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues);
            obj.Addend = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewAddend = ifelse(ParmCodes(end)=='f', obj.Addend,obj.Addend+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewAddend]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = obj.BasisRV.LowerBound + obj.Addend;
            obj.UpperBound = obj.BasisRV.UpperBound + obj.Addend;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.Addend];
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
        
        function thisval = nIthValue(obj,Ith)
            thisval = TransX(obj.BasisRV.nIthValue(Ith));
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Mean(obj.BasisRV) + obj.Addend;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Variance(obj.BasisRV);
        end
        
    end  % methods
    
end  % class AddTrans


