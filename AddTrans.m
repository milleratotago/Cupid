classdef AddTrans < dTransMono
    % AddTrans(BasisRV,Addend): Shifts the BasisRV by the specified additive constant.
    
    properties(SetAccess = protected)
        Addend
    end
    
    methods
        
        function obj=AddTrans(varargin) % BasisDist,Addend
            obj=obj@dTransMono('AddTrans');
            obj.ReviseBounds = false;  % This skips a CDFInverse bounds search and may speed estimation.
            % NEWJEFF: It could also be used in other dTransMono descendants.
            switch nargin
                case 0
                case 2
                    obj.Setup(varargin{:});
                otherwise
                    ME = MException('AddTrans:Constructor','AddTrans constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function Setup(obj,BasisDist,Addend)
            obj.Setup@dTransMono(BasisDist);
            obj.Addend = Addend;
            obj.AddParms(1,'r');
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Addend = newparmvalues(end);
            ResetParms@dTransMono(obj,newparmvalues(1:end-1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            NewAddend = ifelse(ParmCodes(end)=='f', obj.Addend,obj.Addend+0.1);
            obj.ResetParms([obj.BasisRV.ParmValues NewAddend]);
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = obj.Addend;
        end
        
        % Making the next 2 static causes problems:
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = Parms(end);
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = Reals(end);
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = PreTransX + obj.Addend;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = TransX - obj.Addend;
        end
        
        function thisval = PDFScaleFactor(~,~)
            thisval = 1;
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = Mean(obj.BasisRV) + obj.Addend;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = Variance(obj.BasisRV);
        end
        
        function thisval=Skewness(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = Skewness(obj.BasisRV);
        end
        
        function thisval=Kurtosis(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = Kurtosis(obj.BasisRV);
        end
        
        function thisval=CenMoment(obj,I)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = CenMoment(obj.BasisRV,I);
        end
        
    end  % methods
    
end  % class AddTrans


