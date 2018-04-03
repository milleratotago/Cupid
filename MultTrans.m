classdef MultTrans < dTransOf1
    % MultTrans(BasisRV,Multiplier): Rescales the BasisRV by a constant multiplier.
    
    properties(SetAccess = protected)
        Multiplier
        PDFScaler
    end
    
    methods
        
        function obj=MultTrans(varargin)
            obj=obj@dTransOf1('MultTrans');
            obj.NTransParms = 1;
            obj.TransParmCodes = 'r';
            switch nargin
                case 0
                case 2
                    BuildMyBasis(obj,varargin{1});
                    obj.DistType = obj.BasisRV.DistType;
                    obj.NDistParms = obj.BasisRV.NDistParms + 1;
                    obj.DefaultParmCodes = [obj.BasisRV.DefaultParmCodes 'r'];
                    ResetParms(obj,[obj.BasisRV.ParmValues varargin{end}]);
                otherwise
                    ME = MException('MultTrans:Constructor', ...
                        'MultTrans constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues);
            obj.Multiplier = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes);
            obj.BasisRV.PerturbParms(ParmCodes);
            NewMult = ifelse(ParmCodes(end)=='f', obj.Multiplier,obj.Multiplier*1.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewMult]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.TransReverses = obj.Multiplier < 0;
            if obj.Multiplier >= 0
                obj.LowerBound = obj.BasisRV.LowerBound * obj.Multiplier;
                obj.UpperBound = obj.BasisRV.UpperBound * obj.Multiplier;
            else
                obj.LowerBound = obj.BasisRV.UpperBound * obj.Multiplier;
                obj.UpperBound = obj.BasisRV.LowerBound * obj.Multiplier;
            end
            switch obj.DistType
                case 'c'
                    obj.PDFScaler = 1 / abs(obj.Multiplier);
                case 'd'
                    obj.PDFScaler = 1;
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.Multiplier];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX * obj.Multiplier;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX / obj.Multiplier;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
        function thisval = PDFScaleFactor(obj,X)
            thisval = obj.PDFScaler;
        end
        
        function thisval = nIthValue(obj,Ith)
            if obj.Multiplier >= 0
                thisval = TransX(obj.BasisRV.nIthValue(Ith));
            else
                thisval = TransX(obj.BasisRV.nIthValue(obj.NValues-Ith+1));
            end
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


