classdef Poisson < dDiscrete   % NWJEFF: Not vectorized
    % Poisson(mu):  mean equals parameter mu.
    
    % NWJEFF: Use a normal approximation for mu>100; this version gives numerical errors.
    
    properties(SetAccess = protected)
        mu, expmu
        Warned  % Warning given about numerical problems with large mu
    end
    
    methods
        
        function obj=Poisson(varargin)
            obj=obj@dDiscrete('Poisson');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            obj.Warned = false;
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
            ClearBeforeResetParmsD(obj);
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
            if (obj.mu>100)&&(~obj.Warned)
                obj.Warned = true;
                warning('Expect numerical problems with Poisson mu>100; use Normal approximation.');
            end
            obj.expmu = exp(-obj.mu);
            obj.Initialized = true;
            MakeTables(obj);  % also sets lower and upper bounds
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
            thisval = round(X);
        end
        
        function []=MakeTables(obj)
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
            obj.DiscreteX = tempX(1:I);
            obj.DiscretePDF = tempPDF(1:I);
            obj.DiscreteCDF = tempCDF(1:I);
            obj.DiscreteCDF(end) = 1;
            obj.NValues = I;
% NWJEFF: Trim zeros & set NValues?
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
            SetBinEdges(obj);
            obj.StoredTablesInitialized = true;
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



