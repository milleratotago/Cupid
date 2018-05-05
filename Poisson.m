classdef Poisson < dDiscrete   % NWJEFF: Not vectorized
    % Poisson(mu):  mean equals parameter mu.
    
    % NWJEFF: As written, does this _always_ use a stored table?
    
    properties(SetAccess = protected)
        mu, expmu
    end
    
    methods
        
        function obj=Poisson(varargin)
            obj=obj@dDiscrete('Poisson');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('Poisson:Constructor', ...
                        'Poisson constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.StoredTablesInitialized = false;
            obj.mu = newparmvalues;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newmu  = ifelse(ParmCodes(1)=='f', obj.mu, obj.mu*1.1);
            obj.ResetParms(newmu);
        end
        
        function []=ReInit(obj)
            assert(obj.mu>0,'Poisson mu must be > 0.');
            obj.expmu = exp(-obj.mu);
            obj.Initialized = true;
            obj.UseStoredTables = true;
            if obj.UseStoredTables
                MakeTables(obj);  % also sets lower and upper bounds
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.GT2Real(eps,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(eps,Reals(1));
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
%             thisval = zeros(size(X));
%             for i=1:numel(X)
%                 if X(i) < -0.5
%                     thisval(i) = -1;
%                 elseif X(i) < 0.5
%                     thisval(i) = 0;
%                 else
%                     thisval(i) = round(X);
%                 end
                thisval = round(X);
%            end
        end
        
        function thisval=nIthValue(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(I)>0&&max(I)<=obj.NValues,'Requested value at nonexistent position')
            thisval = obj.LowerBound+I-1;
        end
        
        function []=MakeTables(obj)
            obj.LowerBound = 0;
            chunksize = 100;  % grow the arrays in chunks.
            tempX = zeros(1,chunksize);
            tempPDF = zeros(1,chunksize);
            tempCDF = zeros(1,chunksize);
            currentN = chunksize;
            Sum = 0;
            UPower = 1;
            XFac = 1;
            I = 0;
            while Sum<obj.CDFNearlyOne
                if I > currentN
                    tempX = [tempX zeros(1,chunksize)];
                    tempPDF = [tempPDF zeros(1,chunksize)];
                    tempCDF = [tempCDF zeros(1,chunksize)];
                    currentN = currentN + chunksize;
                end
                tempX(I+1) = I;
                tempPDF(I+1) = obj.expmu*UPower/XFac;
                Sum = Sum + tempPDF(I+1);
                tempCDF(I+1) = Sum;
                UPower = UPower * obj.mu;
                I = I+1;
                XFac = XFac * I;
            end
            obj.UpperBound = I;
            obj.NValues = I;
            obj.StoredX = [ -0.01 tempX(1:obj.NValues)];
            obj.StoredPDF = [0 tempPDF(1:obj.NValues)];
            obj.StoredCDF = [0 tempCDF(1:obj.NValues)];
            obj.StoredTablesInitialized = true;
        end

        function thisval=NextValue(obj,X)
            thisval = X + 1;
        end
        
        function thispdf=nPDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            for i=1:numel(X)
                if LegalValue(obj,X(i))
                    thispdf(i) = obj.expmu*obj.mu^X(i)/factorial(X(i));
                end
            end
        end
        
        function thiscdf=nCDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            for i=1:numel(X)
                PSum = 0;
                for I = 0:floor(X(i))
                    PSum = PSum + PDF(obj,I);
                end
                thiscdf(i) = PSum;
            end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.mu;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.mu;
        end
  
        function EstML(obj,Observations,varargin)
            newmu = mean(Observations);
            obj.ResetParms(newmu);
        end
        
    end  % methods
    
end  % class Poisson



