classdef ExGauMn < ExGauss
    % ExGaussian distribution with third parameter = exponential mean (i.e. 1/rate)
    
    properties(SetAccess = protected)
        exmean
    end
    
    methods (Static)

        function parms = MomsToParms(NorMean,NorVar,ExpMean)
            % Return values of 3 distribution parameters yielding specified
            % mean and variance of normal and mean of exponential.
            % Used with ExHelpStartParmsMLE
            parms = [NorMean, sqrt(NorVar), ExpMean];
        end

    end % methods (Static)

    methods
        
        function obj=ExGauMn(varargin)
            obj=obj@ExGauss;
            obj.FamilyName = 'ExGauMn';
            obj.ParmNames{3} = 'exmean';
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExGauMn:Constructor', ...
                        'ExGauMn constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.exmean = newparmvalues(3);
            obj.rate = 1/obj.exmean;
            if ~(obj.sigma>0)
                disp('error');
            end
            if obj.sigma<obj.MinSigma
                error([obj.FamilyName ' sigma is ' num2str(obj.sigma) ' but must be > obj.MinSigma = ' num2str(obj.MinSigma) '.']);
            end
            if obj.exmean<eps
                error(['ExGauMn exmean is ' num2str(obj.exmean) ' but must be > 0.']);
            end
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            if obj.sigma > obj.exmean
                newsigma = 0.95*obj.sigma;
                newmean = 1.05*obj.exmean;
            else
                newsigma = 1.05*obj.sigma;
                newmean = 0.95*obj.exmean;
            end
            newmu    = ifelse(ParmCodes(1)=='f', obj.mu,   1.05*obj.mu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.sigma,newsigma);
            newmean  = ifelse(ParmCodes(3)=='f', obj.exmean, newmean);
            obj.ResetParms([newmu newsigma newmean]);
        end
        
%       Problem: During EstML, ParmsFrom3Mom is called by parent ExGauss & wants real rate
%         function parms = ParmsFrom3Mom(obj,targ1,targ2,targ3)
%             parms = ParmsFrom3Mom@ExGauss(obj,targ1,targ2,targ3);
%             parms(3) = 1 / parms(3);
%         end
        
        function []=MomSetParms(obj,m,s,skew1,ssqr)
            % apply moment-based estimation formulas with this parameterization
            obj.mu = m - s*(skew1/2)^(1/3);
            obj.sigma = sqrt(  ssqr * ( 1 - (skew1/2)^(2/3) )  );
            obj.rate = 1 / (s * (skew1/2)^(1/3));
            obj.exmean = 1 / obj.rate;
        end
        
        function parms = StartParmsMLE(obj,X)
            parms = StartParmsMLE@ExGauss(obj,X);
            parms(3) = 1 / parms(3);
%             HoldParms = obj.ParmValues;
%             parms = ExHelpStartParmsMLE(obj,X);
%             obj.ResetParms(HoldParms);
        end

        function []=EstMom(obj,TargetVals,varargin)
            EstMom@ExGauss(obj,TargetVals,varargin{:});
            obj.exmean = obj.rate;
            obj.rate = 1 / obj.exmean;
            obj.ReInit;
        end

    end
    
end  % class ExGauMn

