classdef Cauchy < dContinuous
    
    properties(SetAccess = protected)
        location, scale,
        OneOverPi, OneOverPiScale
    end
    
    methods
        
        function obj=Cauchy(varargin)
            obj=obj@dContinuous('Cauchy');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.OneOverPi = 1 / pi;
            obj.CDFNearlyZero = 0.001;   % Cut off more than usual of infinite tails
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            obj.IntegralPDFXmuNAbsTol = 100*obj.IntegralPDFXmuNAbsTol;  % For integrating PDF
            obj.IntegralPDFXmuNRelTol = 100*obj.IntegralPDFXmuNRelTol;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Cauchy:Constructor', ...
                        'Cauchy constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.location = newparmvalues(1);
            obj.scale = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newloc   = ifelse(ParmCodes(1)=='f', obj.location, 1.1*obj.location);
            newscale = ifelse(ParmCodes(2)=='f', obj.scale,   0.9*obj.scale);
            obj.ResetParms([newloc newscale]);
        end
        
        function []=ReInit(obj)
            assert(obj.scale>0,'Cauchy scale must be > 0.');
            obj.OneOverPiScale = obj.OneOverPi / obj.scale;
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf=zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            SqrDev = ( (X(InBounds) - obj.location) / obj.scale ).^2;
            thispdf(InBounds) = obj.OneOverPiScale ./ (1 + SqrDev);
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf=zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            ATanDev = atan( (X(InBounds) - obj.location) / obj.scale );
            thiscdf(InBounds) = 0.5 + obj.OneOverPi * ATanDev;
            thiscdf(X>obj.UpperBound) = 1;
        end
        
        function thisval=InverseCDF(obj,P)
            % from Mathematica
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            Ang = pi * (P - 0.5);
            %SinAng = sin(Ang);       % Note: Tan(x) = Sin(x) / Cos(x)
            %CosAng = cos(Ang);
            thisval = obj.location + obj.scale * tan(Ang);
        end
        
        %  The following are theoretically correct but is not true for the bounded approximation being implemented.
        %
        %        function thisval=MGF(obj,Theta)
        %            assert(obj.Initialized,UninitializedError(obj));
        %            thisval = NaN;
        %        end
        %
        %        function thisval=RawMoment(obj,I)
        %            assert(obj.Initialized,UninitializedError(obj));
        %            if I == 0
        %                thisval = 1;
        %            else
        %                thisval = NaN;
        %            end
        %        end
        %
        %        function thisval=CenMoment(obj,I );
        %            assert(obj.Initialized,UninitializedError(obj));
        %            if I == 0
        %                thisval = 1;
        %            else
        %                thisval = NaN;
        %            end
        %        end
        %
        %        function thisval=EstimatesFromMoments(obj,PassMoments,ParmCodes)
        %            ME = MException('Cauchy:EstimatesFromMoments', ...
        %                'Cauchy cannot estimate from moments because they do not exist.');
        %            throw(ME);
        %        end
        
    end  % methods
    
end  % class Cauchy

