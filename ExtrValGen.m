classdef ExtrValGen < dContinuous
    % Generalized Extreme Value distribution with parameters Mu, Sigma>0, Epsilon.
    %   Parameters a.k.a. location, scale, shape
    % Reference: https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
    % Also known as the Fisher-Tippett distribution.
    % Special cases:
    %   ExtrVal1 or Gumbel where Epsilon = 0
    %   ExtrVal2 or Frechet where Epsilon > 0
    %   ExtrVal3 or reversed Weibull where Epsilon < 0
    
    properties(SetAccess = protected)
        Mu, Sigma, Epsilon,
        EpsOverSigma, NegOneOverEps, UsingTypeI, TypeI
        TooSmallAbsEpsilon
    end
    
    methods
        
        function obj=ExtrValGen(varargin)
            obj=obj@dContinuous('ExtrValGen');
            obj.NDistParms = 3;
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.TypeI = ExtrVal1;
            obj.UsingTypeI = false;
            obj.TooSmallAbsEpsilon = sqrt(eps);
            obj.CDFNearlyZero = 1e-8;
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExtrValGen:Constructor', ...
                        'ExtrValGen constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Initialized = false;
            obj.Mu = newparmvalues(1);
            obj.Sigma = newparmvalues(2);
            obj.Epsilon = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newMu = ifelse(ParmCodes(1)=='f', obj.Mu, 1.02*obj.Mu);
            newSigma = ifelse(ParmCodes(2)=='f', obj.Sigma,   0.98*obj.Sigma);
            newEpsilon = ifelse(ParmCodes(3)=='f', obj.Epsilon,   0.98*obj.Epsilon);
            obj.ResetParms([newMu newSigma newEpsilon]);
        end
        
        function ReInit(obj)
            assert(obj.Sigma>0,'ExtrValGen Sigma must be > 0.');
            if abs(obj.Epsilon) <= obj.TooSmallAbsEpsilon
                obj.UsingTypeI = true;
            end
            obj.Initialized = true;
            if obj.UsingTypeI
                ResetParms(obj.TypeI,[obj.Mu obj.Sigma]);
                obj.LowerBound = obj.TypeI.LowerBound;
                obj.UpperBound = obj.TypeI.UpperBound;
            else
                obj.EpsOverSigma = obj.Epsilon / obj.Sigma;
                obj.NegOneOverEps = - 1 / obj.Epsilon;
                if obj.Epsilon > 0
                    obj.LowerBound = obj.Mu - obj.Sigma/obj.Epsilon;
                    % SearchUpper;
                    Increment = 10*obj.Sigma;
                    obj.UpperBound = obj.Mu + Increment;
                    while obj.CDF(obj.UpperBound-1) < obj.CDFNearlyOne
                        obj.UpperBound = obj.Mu + Increment;
                        Increment = 2 * Increment;
                    end
                elseif obj.Epsilon < 0
                    obj.UpperBound =  obj.Mu - obj.Sigma/obj.Epsilon;
                    obj.LowerBound =  obj.Mu + obj.Sigma/obj.Epsilon;
                    % SearchLower;
                    Increment = 10*obj.Sigma;
                    while obj.CDF(obj.LowerBound+1) > obj.CDFNearlyZero
                        obj.LowerBound = obj.Mu - Increment;
                        Increment = 2 * Increment;
                    end
                obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
                obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
                end
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(~,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2)) Parms(3)];
        end
        
        function Parms = RealsToParms(~,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2)) Reals(3)];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if obj.UsingTypeI
                thispdf = PDF(obj.TypeI,X);
            else
                thispdf = zeros(size(X));
                % InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
                OnePlusEpsXmMuOverSigma = 1 + obj.EpsOverSigma * (X(InBounds) - obj.Mu);
                ExpTerm = -OnePlusEpsXmMuOverSigma.^obj.NegOneOverEps;
                thispdf(InBounds) = OnePlusEpsXmMuOverSigma.^(obj.NegOneOverEps-1).*exp(ExpTerm)/obj.Sigma;
                % for i=1:numel(X)
                %     if (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
                %         OnePlusEpsXmMuOverSigma = 1 + obj.EpsOverSigma * (X(i) - obj.Mu);
                %         ExpTerm = -OnePlusEpsXmMuOverSigma^obj.NegOneOverEps;
                %         thispdf(i) = OnePlusEpsXmMuOverSigma^(obj.NegOneOverEps-1)*exp(ExpTerm)/obj.Sigma;
                %     end
                % end
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            if obj.UsingTypeI
                thiscdf = CDF(obj.TypeI,X);
            else
                thiscdf = zeros(size(X));
                thiscdf(X>=obj.UpperBound) = 1;
                OnePlusEpsXmMuOverSigma = 1 + obj.EpsOverSigma * (X(InBounds) - obj.Mu);
                ExpTerm = -OnePlusEpsXmMuOverSigma.^obj.NegOneOverEps;
                thiscdf(InBounds) = exp(ExpTerm);
%                 for i=1:numel(X)
%                     if X(i)<=obj.LowerBound
%                     elseif X(i)>=obj.UpperBound
%                         thiscdf(i) = 1;
%                     else
%                         OnePlusEpsXmMuOverSigma = 1 + obj.EpsOverSigma * (X(i) - obj.Mu);
%                         ExpTerm = -OnePlusEpsXmMuOverSigma^obj.NegOneOverEps;
%                         thiscdf(i) = exp(ExpTerm);
%                     end
%                 end
            end
        end
        
    end  % methods
    
end  % class ExtrValGen

