classdef t < dContinuous
    % t(df)
    
    properties(SetAccess = protected)
        df, LnGamHalfSum, LnGamHalf, SqrtPiDf
        ParentF, NormalApprox
        EvenMomentWarned, NormalConstructed, UsingNormalApprox
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = NumTrans.GT2Real(1,Parms(1));
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = NumTrans.Real2GT(1,Reals(1));
        end
        
    end
    
    methods
        
        function obj=t(varargin)
            obj=obj@dContinuous('t');
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            obj.NDistParms = 1;
            obj.ParentF = F;
            obj.EvenMomentWarned = false;
            obj.NormalConstructed = false;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('t:Constructor', ...
                        't constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.df = VerifyIntegerGE(obj,1,newparmvalues(1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            if ~(ParmCodes(1)=='f')
                obj.ResetParms(obj.df+1);
            end
        end
        
        function []=ReInit(obj)
            obj.UsingNormalApprox = obj.df > 1000;
            if obj.UsingNormalApprox
                if ~obj.NormalConstructed
                    obj.NormalConstructed = true;
                    obj.NormalApprox = Normal(0,1);
                end
                obj.NormalApprox.ResetParms([obj.NormalApprox.Mean obj.NormalApprox.SD]);
                obj.UpperBound = obj.NormalApprox.UpperBound;
                obj.LowerBound = obj.NormalApprox.LowerBound;
            else
                obj.LnGamHalfSum = gammaln( (obj.df + 1)/2 );
                obj.LnGamHalf = gammaln( obj.df / 2 );
                obj.SqrtPiDf = sqrt(pi*obj.df);
                ResetParms(obj.ParentF,[1 obj.df]);
                obj.UpperBound = sqrt(obj.ParentF.UpperBound);
                obj.LowerBound = - obj.UpperBound;
            end
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if obj.UsingNormalApprox
                thispdf = obj.NormalApprox.PDF(X);
            else
                LnOnePlusStuffToPower = log(X(InBounds).^2 / obj.df + 1) * (obj.df + 1) / 2;
                thispdf(InBounds) = exp(obj.LnGamHalfSum - obj.LnGamHalf - LnOnePlusStuffToPower) / obj.SqrtPiDf;
                % for i=1:numel(X)
                %     LnOnePlusStuffToPower = log(X(i).^2 / obj.df + 1) * (obj.df + 1) / 2;
                %     thispdf(i) = exp(obj.LnGamHalfSum - obj.LnGamHalf - LnOnePlusStuffToPower) / obj.SqrtPiDf;
                % end
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            if obj.UsingNormalApprox
                thiscdf(InBounds) = obj.NormalApprox.CDF(X(InBounds));
            else
                XSqr = X.^2;
                thiscdf( InBounds & (X<0) ) = 0.5 - obj.ParentF.CDF(XSqr( InBounds & (X<0) )) / 2;
                thiscdf( InBounds & (X>=0) ) = 0.5 + obj.ParentF.CDF(XSqr( InBounds & (X>=0) )) / 2;
                % for i=1:numel(X)
                %     if X(i) <= obj.LowerBound
                %     elseif X(i) >= obj.UpperBound
                %         thiscdf(i) = 1;
                %     else
                %         XSqr = X(i)^2;
                %         if X(i) < 0
                %             thiscdf(i) = 0.5 - obj.ParentF.CDF(XSqr) / 2;
                %         else
                %             thiscdf(i) = 0.5 + obj.ParentF.CDF(XSqr) / 2;
                %         end
                %     end
                % end
            end
        end
        
        function thisval=Mean(obj)
            thisval = 0;
        end
        
        function thisval=RawSkewness(obj)
            thisval = 0;
        end
        
        % From F: Devroye p 446
        function thisval=Random(obj,varargin)
            if obj.UsingNormalApprox
                thisval = obj.NormalApprox.Random(varargin{:});
            else
                thisval = obj.ParentF.Random(varargin{:}).^0.5;
                u01 = rand(varargin{:});
                s = ones(varargin{:});
                s(u01<.5) = -1;
                thisval = thisval .* s;
            end
        end
        
        function thisval=Variance(obj)    % Mathematica
            if obj.df > 2
                thisval = obj.df / (obj.df - 2);
            else
                thisval = nan;
                if ~obj.EvenMomentWarned
                    obj.EvenMomentWarned = true;
                    warning('2nd & higher even moments do not exist for t distribution with 2 or fewer df.');
                end
            end
        end
        
        function thisval=Kurtosis(obj)    % Mathematica
            if obj.df > 4
                thisval = 6 / (obj.df - 4) + 3;
            else
                thisval = nan;
                if ~obj.EvenMomentWarned
                    obj.EvenMomentWarned = true;
                    warning('4th & higher even moments do not exist for t distribution with 4 or fewer df.');
                end
            end
        end
        
    end  % methods
    
end  % class t






