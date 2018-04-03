classdef PowerTrans < dTransOf1
    % PowerTrans(BasisRV,Power): Power transformation of a BasisRV which must be NONNEGATIVE.
    
    properties(SetAccess = protected)
        Power, InversePower, PowerMinus1
    end
    
    methods
        
        function obj=PowerTrans(varargin)
            obj=obj@dTransOf1('PowerTrans');
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
                    ME = MException('PowerTrans:Constructor', ...
                        'PowerTrans constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues);
            obj.Power = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewPower = ifelse(ParmCodes(end)=='f', obj.Power,obj.Power*1.02);
            obj.ResetParms([obj.BasisRV.ParmValues NewPower]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            assert(obj.BasisRV.LowerBound>=0,'PowerTrans BasisRV must not be negative.');
            % Programming note: To extend this transformation to allow for BasisRVs with negative values, you must at least address the limitations marked *** in this file.
            obj.InversePower = 1 / obj.Power;
            obj.PowerMinus1 = obj.Power - 1;
            % The following 2 lines wouldn't work if obj.BasisRV.obj.LowerBound < 0. ***
            obj.LowerBound = obj.BasisRV.LowerBound^obj.Power;
            obj.UpperBound = obj.BasisRV.UpperBound^obj.Power;
            % Improve bounds to avoid numerical errors
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [obj.Power];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX.^obj.Power;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX.^obj.InversePower;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
        function thisval = PDFScaleFactor(obj,X)
            thisval = ones(size(X));
            switch obj.DistType
                case 'c'
                    PreTransX = X.^obj.InversePower;
                    thisval = thisval ./ abs( obj.Power*PreTransX.^obj.PowerMinus1 );
            end
        end
        
        function thisval = nIthValue(obj,Ith)
            thisval = TransX(obj.BasisRV.nIthValue(Ith));
        end
        
    end  % methods
    
end  % class PowerTrans



