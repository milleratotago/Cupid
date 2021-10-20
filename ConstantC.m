classdef ConstantC < dContinuous
    % Constant defined as a continuous RV
    
    properties(SetAccess = protected)
        value
    end

    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = Parms(1);
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = Reals(1);
        end
        
    end
    
    methods
        
        function obj=ConstantC(varargin)
            obj=obj@dContinuous('ConstantC');
            obj.NDistParms = 1;
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('ConstantC:Constructor', ...
                        'ConstantC constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            % CheckBeforeResetParms(obj,newparmvalues);
            obj.value = newparmvalues(1);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            obj.LowerBound = obj.value;
            obj.UpperBound =  obj.value;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thispdf = zeros(size(X));
            thispdf(X==obj.value) = 1;
        end
        
        function thiscdf=CDF(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thiscdf = zeros(size(X));
            thiscdf(X>=obj.value) = 1;
        end
        
        function thisval=MGF(obj,Theta)
            thisval = exp(Theta * obj.value);
        end
        
        function thisval=InverseCDF(obj,P)
            thisval = ones(size(P)) * obj.value;
        end
        
        function thisval=RawMoment(obj,I)
            thisval = obj.value^I;
        end
        
        function thisval=CenMoment(~,I)
            if I == 0
                thisval = 1;
            else
                thisval = 0;
            end
        end
        
        function thisval=Random(obj,varargin)
            thisval = ones(varargin{:})* obj.value;
        end
        
        function thisval=IntegralXToNxPDF(obj,FromX,ToX,N)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if (FromX <= obj.value) && (obj.value <= ToX)
                thisval = obj.value^N;
            else
                thisval = 0;
            end
        end
        
        function thisval=IntegralX_CToNxPDF(obj,FromX,ToX,C,N)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 0^N;
        end
        
        function thisval=ConditionalRawMoment(obj, FromX, ToX, I)
            if (FromX <= obj.value) && (ToX >= obj.value)
                thisval = RawMoment(obj,I);
            else
                thisval = 0;
            end
        end
        
        function thisval=Hazard(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = zeros(size(X));
            thisval(X==obj.value) = 1;
        end
        
    end  % methods
    
    
end  % class ConstantC


