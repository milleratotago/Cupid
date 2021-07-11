classdef ConstantD < dDiscrete
    % Constant defined as a discrete RV
    
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

            obj.MakeTables;
            SetBinEdges(obj);

            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function []=MakeTables(obj)
            obj.NValues = 1;
            obj.DiscreteX = obj.value;
            obj.DiscretePDF = 1;
            obj.DiscreteCDF = 1;
            obj.StoredTablesInitialized = true;
        end

        function thisval=MGF(obj,Theta)
            thisval = exp(Theta * obj.value);
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
    
    
end  % class ConstantD


