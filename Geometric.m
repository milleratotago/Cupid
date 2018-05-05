classdef Geometric < dDiscrete     % NWJEFF: Not vectorized
    % Geometric(P): Number of trials to first success, including the success (i.e., X=1,2,3,...)
    
    % I think this should always use a table.
    
    properties(SetAccess = protected)
        P, Q, MaxValue
    end
    
    methods
        
        function obj=Geometric(varargin)
            obj=obj@dDiscrete('Geometric');
            obj.DistType = 'd';
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            obj.MaxValue = 1000000;  % Prevents infinite loop in case of P close to 0.
            obj.Grain = 5;  % Need to be very tolerant of CDF machine precision errors.
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('Geometric:Constructor', ...
                        'Geometric constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.P = newparmvalues(1);
            obj.StoredTablesInitialized = false;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newP  = ifelse(ParmCodes(1)=='f', obj.P, obj.P*0.95);
            obj.ResetParms(newP);
        end
        
        function []=ReInit(obj)
            assert(obj.P>0&&obj.P<1,'Geometric P must be 0-1.');
            obj.Q = 1 - obj.P;
            obj.LowerBound = 1;
            ThisI = 1;
            ThisPDF = obj.P;
            while (ThisPDF>obj.CDFNearlyZero) && (ThisI<obj.MaxValue)
                ThisI = ThisI + 1;
                ThisPDF = ThisPDF*obj.Q;
            end
            obj.UpperBound = ThisI;
            obj.NValues = round(obj.UpperBound - obj.LowerBound) + 1;
            obj.Initialized = true;
            if obj.UseStoredTables
                obj.MakeTables;
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.Bounded2Real(0,1,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2Bounded(0,1,Reals(1));
        end
        
        function thisval=LegalValue(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = zeros(size(X));
            for i=1:numel(X)
                if (abs(round(X(i))-X(i))<eps) && (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
                    thisval(i) = true;
                end
            end
        end
        
        function thisval=NearestLegal(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = zeros(size(X));
            for i=1:numel(X)
                if X(i) < 1.5
                    thisval(i) = 1;
                elseif X(i)>obj.UpperBound-0.5
                    thisval(i) = obj.UpperBound;
                else
                    thisval(i) = round(X(i));
                end
            end
        end
        
        function thisval=nIthValue(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(I)>=obj.LowerBound&&max(I)<=obj.UpperBound,'Requested value at nonexistent position')
            thisval = obj.LowerBound+I-1;
        end
        
        function thisval=IthValue(obj,I)  % NEWJEFF: Does not handle array of I's
            thisval = I - 1 + obj.LowerBound;
        end
        
        function thisval=NextValue(obj,X)
            thisval = X + 1;
        end
        
        function thispdf=nPDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            Legal = obj.LegalValue(X);
            Xi = round(X);
            InBounds = (Xi>=obj.LowerBound) & (Xi<=obj.UpperBound);
            thispdf(Legal & InBounds) = obj.Q.^(Xi(Legal & InBounds)-1) * obj.P;
        end
        
        function thiscdf=nCDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(InBounds) = 1 - (1-obj.P).^floor(X(InBounds));
            thiscdf(X>obj.UpperBound) = 1;
            % zeros(size(X));
            % for i=1:numel(X)
            %     PSum = 0;
            %     for I = 0:floor(X(i))
            %         PSum = PSum + PDF(obj,I);
            %     end
            %     thiscdf(i) = PSum;
            % end
        end% 
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1/obj.P;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (1-obj.P)/obj.P^2;
        end
  
        function thisval=Skewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (2-obj.P)/sqrt(1-obj.P);
        end
  
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 9 + obj.P^2/(1-obj.P);
        end
  
        function s=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            meanObs = mean(Observations);
            ResetParms(obj,1/meanObs);
            BuildMyName(obj);
            s=obj.StringName;
        end
        
        function thisval=Random(obj,varargin)
            thisval = geornd(obj.P,varargin{:}) + 1;  % geornd is number of failures before first success
        end
        
    end  % methods
    
end  % class Geometric



