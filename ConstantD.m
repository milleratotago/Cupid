classdef ConstantD < dDiscrete
    % Constant defined as a discrete RV
    
    properties(SetAccess = protected)
        value
    end
    
    methods
        
        function obj=ConstantD(varargin)
            obj=obj@dDiscrete('ConstantD');
            obj.NDistParms = 1;
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('ConstantD:Constructor', ...
                        'ConstantD constructor needs 0 or 1 arguments.');
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

            % MakeTables:
            obj.NValues = 1;
            obj.DiscreteX = obj.value;
            obj.DiscretePDF = 1;
            obj.DiscreteCDF = 1;

            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = Parms(1);
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = Reals(1);
        end
        
        function thisval=MGF(obj,Theta)
            thisval = exp(Theta * obj.value);
        end
        
        function thisval=RawMoment(obj,I)
            thisval = obj.value^I;
        end
        
        function thisval=CenMoment(obj,I)
            if I == 0
                thisval = 1;
            else
                thisval = 0;
            end
        end
        
        function thisval=Random(obj,varargin)
            thisval = ones(varargin{:})* obj.value;
        end
        
        function thisval=ConditionalRawMoment(obj, FromX, ToX, I)
            if (FromX <= obj.value) && (ToX >= obj.value)
                thisval = RawMoment(obj,I);
            else
                thisval = 0;
            end
        end
        
    end  % methods
    
    
end  % class ConstantD


