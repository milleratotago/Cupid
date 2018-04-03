classdef TruncatedXhi < TruncatedX
    % TruncatedXhi(BasisRV,UpperCutoffX) produces the truncated version of the BasisRV,
    % truncating at the high end with the indicated X cutoff.
    
    methods
        
        function obj=TruncatedXhi(varargin)
            obj=obj@TruncatedX;
            obj.ThisFamilyName = 'TruncatedXhi';
            obj.NTransParms = 1;
            obj.TransParmCodes = 'r';
            switch nargin
                case 0
                case 2
                    BuildMyBasis(obj,varargin{1});
                    ResetParms(obj,[ParmValues(obj.BasisRV) varargin{end}]);
                otherwise
                    ME = MException('TruncatedXhi:Constructor', ...
                        'TruncatedXhi constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues(1:end-1));
            obj.UpperCutoffX = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewUpper = ifelse(ParmCodes(end)=='f', obj.UpperCutoffX,obj.UpperCutoffX+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewUpper]);
        end
        
        function []=ReInit(obj)
            obj.LowerCutoffX = obj.BasisRV.LowerBound;
            obj.LowerBound = obj.BasisRV.LowerBound;
            if obj.UpperCutoffX <= obj.BasisRV.UpperBound
                obj.UpperBound = obj.UpperCutoffX;
            else
                obj.UpperBound = obj.BasisRV.UpperBound;
            end
            obj.NDistParms = obj.BasisRV.NDistParms + 1;
            obj.DefaultParmCodes = [obj.BasisRV.DefaultParmCodes 'f'];
            obj.DistType = obj.BasisRV.DistType;
            obj.NValues = obj.BasisRV.NValues;
            obj.IthValOffset = 0;
            switch obj.DistType
                case 'd'
                    obj.XTolerance = obj.XNearlyZero;
                    while IthValue(obj.BasisRV,obj.NValues) > obj.UpperBound
                        obj.NValues = obj.NValues - 1;
                    end
                    StillLooking = true;
                    while StillLooking
                        obj.IthValOffset = obj.IthValOffset + 1;
                        StillLooking = IthValue(obj.BasisRV,obj.IthValOffset) < obj.LowerBound;
                    end
                    obj.IthValOffset = obj.IthValOffset - 1;
                    obj.NValues = obj.NValues - obj.IthValOffset;
                case 'c'
                    obj.XTolerance = 0;
            end
            obj.LowerCutoffP = 0;
            obj.UpperCutoffP = CDF(obj.BasisRV,obj.UpperBound+obj.XTolerance);
            obj.UnconditionalP = obj.UpperCutoffP - obj.LowerCutoffP;
            assert(obj.UnconditionalP>0,'TruncatedXhi distribution must include probability > 0.');
            
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = ParmValues(obj,varargin)
            parmvals = [obj.BasisRV.ParmValues obj.UpperCutoffX];
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.UpperCutoffX;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
    end  % methods
    
end  % class TruncatedXhi



