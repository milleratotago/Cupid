classdef ExWaldMn < ExWald
    % ExWaldMn(mu, sigma, barrierA, expmean): Like ExWald except last parameter is mean of exponential instead of rate.
    
    properties(SetAccess = protected)
        expmean
    end
    
    methods
        
        function obj=ExWaldMn(varargin)
            obj=obj@ExWald;
            obj.FamilyName = 'ExWaldMn';
            obj.ParmNames{4} = 'expmean';
            switch nargin
                case 0
                case 4
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExWaldMn:Constructor', ...
                        'ExWaldMn constructor needs 0 or 4 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.barrierA = newparmvalues(3);
            obj.expmean = newparmvalues(4);
            if obj.mu<=0
                error('ExWaldMn mu must be > 0.');
            end
            if obj.sigma<=0
                error('ExWaldMn sigma must be > 0.');
            end
            if obj.barrierA<=0
                error('ExWaldMn barrierA must be > 0.');
            end
            if obj.expmean<=0
                error('ExWaldMn exponential mean must be > 0.');
            end
            obj.rate = 1/obj.expmean;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newmu      = ifelse(ParmCodes(1)=='f', obj.mu, 1.05 * obj.mu);
            newsigma   = ifelse(ParmCodes(2)=='f', obj.sigma,1.05 * obj.sigma);
            newbarrier = ifelse(ParmCodes(3)=='f', obj.barrierA,1.05 * obj.barrierA);
            newexpmean    = ifelse(ParmCodes(4)=='f', obj.expmean,1.05*obj.expmean);
            obj.ResetParms([newmu newsigma newbarrier newexpmean]);
        end
        
    end  % methods
    
end  % class ExWaldMn

