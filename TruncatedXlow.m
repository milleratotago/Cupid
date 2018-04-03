classdef TruncatedXlow < TruncatedX
    % TruncatedXlow(BasisRV,LowerCutoffX) produces the truncated version of the BasisRV,
    % truncating at the low end with the indicated X cutoff.
    
    methods
        
        function obj=TruncatedXlow(varargin)
            obj=obj@TruncatedX;
            obj.ThisFamilyName = 'TruncatedXlow';
            obj.NTransParms = 1;
            obj.TransParmCodes = 'r';
            switch nargin
                case 0
                case 2
                    BuildMyBasis(obj,varargin{1});
                    ResetParms(obj,[ParmValues(obj.BasisRV) varargin{end}]);
                otherwise
                    ME = MException('TruncatedXlow:Constructor', ...
                        'TruncatedXlow constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues(1:end-1));
            obj.LowerCutoffX = newparmvalues(end);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewLower = ifelse(ParmCodes(end)=='f', obj.LowerCutoffX,obj.LowerCutoffX-0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewLower]);
        end
        
        function []=ReInit(obj)
            if obj.LowerCutoffX <= obj.BasisRV.LowerBound
                obj.LowerBound = obj.BasisRV.LowerBound;
            else
                obj.LowerBound = obj.LowerCutoffX;
            end
            obj.UpperCutoffX = obj.BasisRV.UpperBound;
            obj.UpperBound = obj.BasisRV.UpperBound;
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
            obj.LowerCutoffP = CDF(obj.BasisRV,obj.LowerBound-obj.XTolerance);
            obj.UpperCutoffP = 1;
            obj.UnconditionalP = obj.UpperCutoffP - obj.LowerCutoffP;
            assert(obj.UnconditionalP>0,'TruncatedXlow distribution must include probability > 0.');
            
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = ParmValues(obj,varargin)
            parmvals = [obj.BasisRV.ParmValues obj.LowerCutoffX];
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.LowerCutoffX;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
    end  % methods
    
end  % class TruncatedXlow



