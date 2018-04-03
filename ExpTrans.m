classdef ExpTrans < dTransOf1
    % ExpTrans(BasisRV): Exponential transformation of a BasisRV.
    
    properties(SetAccess = protected)
    end
    
    methods
        
        function obj=ExpTrans(varargin)
            obj=obj@dTransOf1('ExpTrans');
            obj.NTransParms = 0;
            obj.TransParmCodes = '';
            switch nargin
                case 0
                case 1
                    BuildMyBasis(obj,varargin{1});
                    obj.DistType = obj.BasisRV.DistType;
                    obj.NDistParms = obj.BasisRV.NDistParms;
                    obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
                    ResetParms(obj,obj.BasisRV.ParmValues);
                otherwise
                    ME = MException('ExpTrans:Constructor', ...
                        'ExpTrans constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            obj.ResetParms(obj.BasisRV.ParmValues);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = PreTransToTrans(obj,obj.BasisRV.LowerBound);
            obj.UpperBound = PreTransToTrans(obj,obj.BasisRV.UpperBound);
            % Improve bounds to avoid numerical errors
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = exp(PreTransX);
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = log(TransX);
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = [];
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = [];
        end
        
        function thisval = PDFScaleFactor(obj,X)
            thisval=ones(size(X));
            switch obj.DistType
                case 'c'
                    thisval = thisval ./ X;
            end
        end
        
        function thisval = nIthValue(obj,Ith)
            thisval = TransX(obj.BasisRV.nIthValue(Ith));
        end
        
    end  % methods
    
end  % class ExpTrans

