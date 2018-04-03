classdef LinearTrans < dTransOf1
    % LinearTrans(BasisRV,Multiplier,Addend) produces the linear transformation BasisRV*Multiplier+Addend.
    
    properties(SetAccess = protected)
        Multiplier, Addend
        PDFScaler
    end
    
    methods
        
        function obj=LinearTrans(varargin)
            obj=obj@dTransOf1('LinearTrans');
            obj.NTransParms = 2;
            obj.TransParmCodes = 'rr';
            switch nargin
                case 0
                case 3
                    BuildMyBasis(obj,varargin{1});
                    obj.DistType = obj.BasisRV.DistType;
                    obj.NDistParms = obj.BasisRV.NDistParms + 2;
                    obj.DefaultParmCodes = [obj.BasisRV.DefaultParmCodes 'rr'];
                    ResetParms(obj,[obj.BasisRV.ParmValues varargin{end-1} varargin{end}]);
                otherwise
                    ME = MException('LinearTrans:Constructor', ...
                        'LinearTrans constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues);
            obj.Multiplier = newparmvalues(end-1);
            obj.Addend = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewMultiplier = ifelse(ParmCodes(end-1)=='f', obj.Multiplier,0.9*obj.Multiplier);
            NewAddend = ifelse(ParmCodes(end)=='f', obj.Addend,obj.Addend+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewMultiplier NewAddend]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.TransReverses = obj.Multiplier < 0;
            if obj.Multiplier >= 0
                obj.LowerBound = obj.BasisRV.LowerBound * obj.Multiplier + obj.Addend;
                obj.UpperBound = obj.BasisRV.UpperBound * obj.Multiplier + obj.Addend;
            else
                obj.LowerBound = obj.BasisRV.UpperBound * obj.Multiplier + obj.Addend;
                obj.UpperBound = obj.BasisRV.LowerBound * obj.Multiplier + obj.Addend;
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
            parmvals = [obj.Multiplier obj.Addend];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX * obj.Multiplier + obj.Addend;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = (TransX - obj.Addend) / obj.Multiplier;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end-1:end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end-1:end);
        end
        
        function thisval = PDFScaleFactor(obj,~)
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
            thisval = Mean(obj.BasisRV) * obj.Multiplier + obj.Addend;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Variance(obj.BasisRV) * obj.Multiplier^2;
        end
        
    end  % methods
    
end  % class LinearTrans


