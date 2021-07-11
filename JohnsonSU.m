classdef JohnsonSU < dContinuous
    % Johnson SU distribution with parameters Location,Scale,Alpha1>0,Alpha2>0
    % https://en.wikipedia.org/wiki/Johnson's_SU-distribution calls these four parameters (in order):
    %    Xi (squiggly E), lambda (approx -\), gamma (squiggly y), delta (dangly o)
    % Warning: Parameter estimation does not work well (now skipped in utJohnsonSU).
    
    properties(SetAccess = protected)
        Location, Scale, Alpha1, Alpha2
        PDFCon, MyZExtreme, Standard_Normal
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2)) Parms(3) NumTrans.GT2Real(eps,Parms(4))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2)) Reals(3) NumTrans.Real2GT(eps,Reals(4))];
        end
        
    end
    
    methods
        
        function obj=JohnsonSU(varargin)
            obj=obj@dContinuous('JohnsonSU');
            obj.NDistParms = 4;
            obj.ParmTypes = 'rrrr';
            obj.DefaultParmCodes = 'rrrr';
            obj.Standard_Normal = Normal(0,1);
            obj.MyZExtreme = 6;  % May need adjustment depending on Alpha2
            switch nargin
                case 0
                case 4
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('JohnsonSU:Constructor', ...
                        'JohnsonSU constructor needs 0 or 4 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Location = newparmvalues(1);
            obj.Scale = newparmvalues(2);
            obj.Alpha1 = newparmvalues(3);
            obj.Alpha2 = newparmvalues(4);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newlocation = ifelse(ParmCodes(1)=='f', obj.Location,obj.Location+.01);
            newscale  = ifelse(ParmCodes(2)=='f', obj.Scale, 1.01*obj.Scale);
            newalpha1 = ifelse(ParmCodes(3)=='f', obj.Alpha1, 0.99*obj.Alpha1);
            newalpha2 = ifelse(ParmCodes(4)=='f', obj.Alpha2, 1.01*obj.Alpha2);
            obj.ResetParms([newlocation newscale newalpha1 newalpha2]);
        end
        
        function []=ReInit(obj)
            assert(obj.Scale>0,'JohnsonSU Scale must be > 0.');
            assert(obj.Alpha2>0,'JohnsonSU Alpha2 must be > 0.');
            obj.PDFCon = obj.Alpha2/(obj.Scale*sqrt(2*pi));
            obj.LowerBound = obj.Location + obj.Scale * sinh( (-obj.MyZExtreme-obj.Alpha1)/obj.Alpha2 );
            obj.UpperBound = obj.Location + obj.Scale * sinh( ( obj.MyZExtreme-obj.Alpha1)/obj.Alpha2 );
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.PDFCon * 1 ./ sqrt(1 + ((X(InBounds)-obj.Location)/obj.Scale).^2) ...
                .* exp(-0.5*(obj.Alpha1+obj.Alpha2*asinh((X(InBounds)-obj.Location)/obj.Scale)).^2);
            %            for i=1:numel(X)
            %                if (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
            %                    X(i) = (X(i) - obj.Location) / obj.Scale;
            %                    SqrtXSqrPlus1 = sqrt( X(i)^2+1 );
            %                    Dif = X(i) + SqrtXSqrPlus1;
            %                    if Dif <= 0
            %                        Dif = obj.XNearlyZero;  % Numerical errors sometimes give Dif<0. :(
            %                    end
            %                    Z = obj.Alpha1 + obj.Alpha2 * log(Dif);
            %                    NorPDF = PDF(obj.Standard_Normal,Z);
            %                    thispdf(i) = obj.Alpha2 / (obj.Scale*SqrtXSqrPlus1) * NorPDF;
            %                 end
            %             end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = obj.Standard_Normal.CDF( obj.Alpha1 + obj.Alpha2*asinh( (X(InBounds)-obj.Location)/obj.Scale ) );   % asinh = inverse of sinh
            %           for i=1:numel(X)
            %               if X(i)<=obj.LowerBound
            %               elseif X(i)>=obj.UpperBound
            %                   thiscdf(i) = 1;
            %               else
            %                   SqrtXSqrPlus1 = sqrt( X(i)^2+1 );
            %                   Dif = X(i) + SqrtXSqrPlus1;
            %                   if Dif <= 0
            %                       Dif = obj.XNearlyZero;  % Numerical errors sometimes give Dif<0. :(
            %                   end
            %                   Z = obj.Alpha1 + obj.Alpha2 * log( Dif );
            %                   thiscdf(i) = CDF(obj.Standard_Normal,Z);
            %               end
            %           end
        end
        
        function thisval=Mean(obj)  % Wikipedia
            thisval = obj.Location - obj.Scale * exp(obj.Alpha2^(-2)/2) * sinh(obj.Alpha1/obj.Alpha2);
        end
        
        function thisval=Variance(obj)  % Wikipedia
            thisval = obj.Scale^2 / 2 * (exp(obj.Alpha2^(-2)) - 1) * (exp(obj.Alpha2^(-2)) * cosh(2*obj.Alpha1/obj.Alpha2) + 1);
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            Z = Random(obj.Standard_Normal,varargin{:});
            thisval = obj.Location + obj.Scale .* sinh( (Z-obj.Alpha1)/obj.Alpha2 );
        end
        
    end  % methods
    
end  % class JohnsonSU


