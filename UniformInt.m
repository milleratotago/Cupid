classdef UniformInt < dDiscrete
    % UniformInt(Low,High): equally likely integers from Low to High
    
    properties(SetAccess = protected)
        LowInt, HighInt, PrValue
    end
    
    methods
        
        function obj=UniformInt(varargin)
            obj=obj@dDiscrete('UniformInt');
            obj.DistType = 'd';
            obj.NDistParms = 2;
            obj.ParmTypes = 'ii';
            obj.DefaultParmCodes = 'ii';
            obj.UseStoredTables = false;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('UniformInt:Constructor', ...
                        'UniformInt constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.StoredTablesInitialized = false;
            obj.LowInt = floor(newparmvalues(1));
            obj.HighInt = floor(newparmvalues(2));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newLo  = ifelse(ParmCodes(1)=='f', obj.LowInt, obj.LowInt-1);
            newHi  = ifelse(ParmCodes(2)=='f', obj.HighInt, obj.HighInt+1);
            obj.ResetParms([newLo newHi]);
        end
        
        function []=ReInit(obj)
            assert(obj.LowInt<obj.HighInt,'UniformInt requires Low<High.');
            obj.LowerBound = obj.LowInt;
            obj.UpperBound = obj.HighInt;
            obj.NValues = obj.UpperBound - obj.LowerBound + 1;
            obj.PrValue = 1 / obj.NValues;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % It is tempting to enforce the restriction that Parms(2) has to be greater than Parms(1), like this:
            % Reals = [Parms(1) GTToAnyReal(Parms(1),Parms(2))];
            % Unfortunately, this creates problems if you want to estimate Parms(1) while holding Parms(2) fixed.
            % It may be better to use a different parameterization, e.g. Parms(1) = bottom & Parms(2) = width,
            % so that these parameters are independent.
            Reals = Parms;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            % Parms = [Reals(1) AnyRealToGT(Reals(1),Reals(2))];
            Parms = Reals;
        end
        
        function thisval=LegalValue(obj,X)
            thisval = zeros(size(X));
            for i=1:numel(X)
                if (floor(X(i))==ceil(X(i))) && (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
                    thisval(i) = true;
                end
            end
        end
        
        function thisval=NearestLegal(obj,X)
            thisval = round(X);
        end
        
        function thispdf=nPDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = ones(size(X))*obj.PrValue;
            Illegal = ~obj.LegalValue(X);
            thispdf(Illegal) = 0;
        end
        
        function thiscdf=nCDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            thiscdf(X>=obj.LowerBound) = (floor(X(X>=obj.LowerBound))-obj.LowerBound + 1) / obj.NValues;
            thiscdf(thiscdf>1) = 1;
            % for i=1:numel(X)
            %     if (X(i) < obj.LowerBound)
            %     elseif (X(i) < obj.UpperBound)
            %         thiscdf(i) = (floor(X(i))-obj.LowerBound + 1)/obj.NValues;
            %     else
            %         thiscdf(i) = 1;
            %     end
            % end
        end
        
        function thisval=InverseCDF(obj,P)
            thisval = ceil(P*obj.NValues-obj.Grain*eps(P*obj.NValues)) + obj.LowerBound - 1;
        end
        
        function thisval=nIthValue(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(I)>0&&max(I)<=obj.NValues,'Requested value at nonexistent position')
            thisval = I - 1 + obj.LowerBound;
        end
        
        function thisval=NextValue(obj,X)
            thisval = X + 1;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (obj.LowerBound + obj.UpperBound) / 2;
        end
        
        function []=EstimatesFromObservations(obj,PassObservations,ParmCodes)
            if ParmCodes(1) == 'f'
                ObsMin = obj.LowerBound;
            else
                ObsMin = min(PassObservations);
            end
            if ParmCodes(2) == 'f'
                ObsMax = obj.UpperBound;
            else
                ObsMax = max(PassObservations);
            end
            Init(obj,ObsMin,ObsMax);
        end
        
        function thisval=Random(obj,varargin)
            thisval = floor(rand(varargin{:})/obj.PrValue) + obj.LowerBound;
        end
        
    end  % methods
    
end  % class UniformInt



