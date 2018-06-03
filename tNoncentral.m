classdef tNoncentral < dContinuous
    % tNoncentral(df,noncentrality)
    % For a one-sample t-test, the noncentrality parameter is
    %   $\sqrt{n}\theta$, where $\theta$ is the true mean/sigma.
    % For a two-sample t-test with sample sizes n1 and n2, the noncentrality parameter is
    %   $\sqrt{n1*n2/(n1+n2)}\theta$, where $\theta$ is the true (mean1-mean2)/sigma.
    
    % This version uses built-in MATLAB functions nctpdf & nctcdf.
    % In storage is a saved older version converted from Pascal, but it had problems.
    
    properties(SetAccess = protected)
        df, noncen
    end
    
    methods
        
        function obj=tNoncentral(varargin)
            obj=obj@dContinuous('tNoncentral');
            obj.ParmTypes = 'ir';
            obj.DefaultParmCodes = 'fr';
            obj.NDistParms = 2;
            obj.CDFNearlyZero = 1e-8;
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('tNoncentral:Constructor', ...
                        'tNoncentral constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.df = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.noncen = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newdf   = ifelse(ParmCodes(1)=='f', obj.df, obj.df+1);
            newnoncen = ifelse(ParmCodes(2)=='f', obj.noncen,1.1*obj.noncen);
            obj.ResetParms([newdf newnoncen]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = -1000;
            obj.UpperBound = 1000;
            
            % Search for lower bound:
            TryValue = obj.noncen;
            LastCDF = 2;
            ThisCDF = 1;
            while (ThisCDF > obj.CDFNearlyZero) && (LastCDF-ThisCDF > eps) && (TryValue > -Inf)
                LastCDF = ThisCDF;
                TryValue = TryValue - 1;
                ThisCDF = obj.CDF(TryValue);
            end
            obj.LowerBound = TryValue;
            if ThisCDF > obj.CDFNearlyZero
                warning(['tNoncentral CDF lower asymptote >0 with df, noncen = ' num2str(obj.df) ',' num2str(obj.noncen)])
            end
            
            % Search for upper bound:
            TryValue = obj.noncen;
            ThisCDF = -0.5;
            LastCDF = -1.5;
            while (ThisCDF < obj.CDFNearlyOne) && (ThisCDF-LastCDF > eps) && (TryValue < Inf)
                LastCDF = ThisCDF;
                TryValue = TryValue + 1;
                ThisCDF = obj.CDF(TryValue);
            end
            obj.UpperBound = TryValue;
            if ThisCDF < obj.CDFNearlyOne
                warning(['tNoncentral CDF asymptote ' num2str(ThisCDF) ' at X= ' num2str(obj.UpperBound) ' with df, noncen = '  num2str(obj.df) ' ' num2str(obj.noncen)])
            end
            
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) Parms(2)];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) Reals(2)];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = nctpdf(X(InBounds),obj.df,obj.noncen);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = nctcdf(X(InBounds),obj.df,obj.noncen);
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.noncen * sqrt(obj.df/2) * gamma((obj.df-1)/2) / gamma(obj.df/2);
        end
        
        function thisval=Variance(obj)    % Mathematica
            if obj.df <= 2
                thisval = nan;
            else
                thisval = obj.df*(1+obj.noncen^2)/(obj.df-2) - obj.noncen^2*obj.df/2*(  gamma((obj.df-1)/2) / gamma(obj.df/2))^2;
            end
        end
        
    end  % methods
    
end  % class tNoncentral



