classdef ExGauRatio < ExGauss
    % ExGaussian distribution with third parameter = ratio of exponential mn&sd to normal sd
    
    properties(SetAccess = protected)
        ratio
    end
    
    methods
        
        function obj=ExGauRatio(varargin)   % Constructor
            obj=obj@ExGauss;
            obj.FamilyName = 'ExGauRatio';
            obj.ParmNames{3} = 'ratio';
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExGauRatio:Constructor', ...
                        'ExGauRatio constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.ratio = newparmvalues(3);
            obj.rate = 1 / (obj.ratio*obj.sigma);
            if obj.sigma<=0
                error('ExGauRatio sigma must be > 0.');
            end
            if obj.ratio<=0
                error('ExGauRatio ratio must be > 0.');
            end
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            if obj.sigma > obj.ratio*obj.mu
                newsigma = 0.95*obj.sigma;
                newratio = 1.05*obj.ratio;
            else
                newsigma = 1.05*obj.sigma;
                newratio = 0.95*obj.ratio;
            end
            newmu    = ifelse(ParmCodes(1)=='f', obj.mu,    1.05*obj.mu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.sigma, newsigma);
            newratio = ifelse(ParmCodes(3)=='f', obj.ratio, newratio);
            obj.ResetParms([newmu newsigma newratio]);
        end
        
        function []=MomSetParms(obj,m,s,skew1,ssqr)
            % apply moment-based estimation formulas with this parameterization
            obj.mu = m - s*(skew1/2)^(1/3);
            obj.sigma = sqrt(  ssqr * ( 1 - (skew1/2)^(2/3) )  );
            obj.rate = 1 / (s * (skew1/2)^(1/3));
            obj.ratio = (1 / obj.rate)/obj.sigma;
        end

%       Problem: During EstML, ParmsFrom3Mom is called by parent ExGauss & wants real rate
%         function parms = ParmsFrom3Mom(obj,targ1,targ2,targ3)
%             parms = ParmsFrom3Mom@ExGauss(obj,targ1,targ2,targ3);
%             parms(3) = 1 / (parms(2) * parms(3));
%         end
        
    end
    
end  % class ExGauRatio

