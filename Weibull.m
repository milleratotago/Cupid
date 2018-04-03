classdef Weibull < dContinuous
    % Weibull(scale>0,power>0,origin)
    
    % It would be nice to have a way to put an upper limit on the OriginParm
    % for estimation problems to make sure that all known observations exist
    % in the distribution (e.g., MLE & percentile-based estimation).
    % More could be be done for OrderIID statistics from Weibull (see p. 254-5).
    %  The best way to do this might to treat it as a new distribution: WeibulO.
    
    properties(SetAccess = protected)
        scale, power, origin, invpower, Standard_Exponential
    end
    
    methods
        
        function obj=Weibull(varargin)
            obj=obj@dContinuous('Weibull');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.CDFNearlyZero = 0.1e-7;  % Numerical problems in tails
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            obj.Standard_Exponential = Exponential(1);
            %           It may be appropriate to relax these for very small or large values of the Power parameter.
            %            obj.IntegralPDFXNAbsTol = 10*obj.IntegralPDFXNAbsTol;  % For integrating PDF
            %            obj.IntegralPDFXNRelTol = 10*obj.IntegralPDFXNRelTol;
            %            obj.IntegralPDFXmuNAbsTol = 10*obj.IntegralPDFXmuNAbsTol;  % For integrating PDF
            %            obj.IntegralPDFXmuNRelTol = 10*obj.IntegralPDFXmuNRelTol;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Weibull:Constructor', ...
                        'Weibull constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.scale = newparmvalues(1);
            obj.power = newparmvalues(2);
            obj.origin = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Careful: the perturbed Origin must be smaller or the data points will have 0 pdf values!
            newscale  = ifelse(ParmCodes(1)=='f', obj.scale, 1.05*obj.scale);
            newpower  = ifelse(ParmCodes(2)=='f', obj.power, 0.95*obj.power);
            neworigin = ifelse(ParmCodes(3)=='f', obj.origin, obj.origin+.1);
            obj.ResetParms([newscale newpower neworigin]);
        end
        
        function []=ReInit(obj)
            assert(obj.scale>0,'Weibull scale must be > 0.');
            assert(obj.power>0,'Weibull power must be > 0.');
            obj.invpower = 1 / obj.power;
            obj.Initialized = true;
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2)) Parms(3)];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) Reals(3)];
        end
        
        function thispdf=PDF(obj,X)  % Johnson & Kotz 1970, p. 250
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            Y = (X(InBounds) - obj.origin) / obj.scale;
            Y2PwrM1 = Y.^(obj.power-1);
            thispdf(InBounds) = obj.power / obj.scale * Y2PwrM1 .* exp( -(Y2PwrM1.*Y) );
            % for iel=1:numel(X)
            %     if InBounds(iel)
            %         Y = (X(iel) - obj.origin) / obj.scale;
            %         Y2PwrM1 = Y^(obj.power-1);
            %         thispdf(iel) = obj.power / obj.scale * Y2PwrM1 * exp( -(Y2PwrM1*Y) );
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)  % Johnson & Kotz 1970, p. 252
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            Y = (X(InBounds) - obj.origin) / obj.scale;
            Y2Pwr = Y.^obj.power;
            thiscdf(InBounds) = 1 - exp( -Y2Pwr );
            % for iel=1:numel(X)
            %     if InBounds(iel)
            %         Y = (X(iel) - obj.origin) / obj.scale;
            %         Y2Pwr = Y^obj.power;
            %         thiscdf(iel) = 1 - exp( -Y2Pwr );
            %     end
            % end
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            InvY = -log(1 - P(InBounds));
            InvY = InvY.^obj.invpower;
            thisval(InBounds) = InvY * obj.scale + obj.origin;
        end
        
        function thisval=Random(obj,varargin)
            % Devroye states that the standard Weibull is the density of
            %  E^(1/a) where E is an exponential random variable.
            %  This also checks with p. 250 of Johnson & Kotz.
            thisval = Random(obj.Standard_Exponential,varargin{:}).^obj.invpower * obj.scale + obj.origin;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = gamma(1 + 1 / obj.power) * obj.scale + obj.origin;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            Gam1 = gamma(1 + 1 / obj.power);
            Gam2 = gamma(1 + 2 / obj.power);
            thisval = (Gam2 - Gam1^2) * obj.scale^2;
        end
        
    end  % methods
    
end  % class Weibull

