classdef RNGammaMS < RNGamma
    % RNGammaMS(mu,sigma) is an RNGamma specified by its mean and sd.

    % WARNING: Check for numerical problems with strange mu,sigma targets.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, sigma
    end
    
    
    methods
        
        function obj=RNGammaMS(varargin)
            obj=obj@RNGamma;
            obj.ThisFamilyName = 'RNGammaMS';
            obj.ParmNames{1} = 'mu';
            obj.ParmNames{2} = 'sigma';
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('RNGammaMS:Constructor', ...
                        'RNGammaMS constructor must receive 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            assert(obj.mu>0,'RNGammaMS mu must be > 0.');
            assert(obj.sigma>0,'RNGammaMS sigma must be > 0.');
            tempN = (obj.mu/obj.sigma)^2;
            tempRate = tempN / obj.mu;
            HoldNameBuilding = obj.NameBuilding;
            obj.NameBuilding = false;
            ResetParms@RNGamma(obj,[tempN tempRate]);
            obj.NameBuilding = HoldNameBuilding;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            NewMean = ifelse(ParmCodes(1)=='f',obj.mu,1.05*obj.mu);
            NewSigma = ifelse(ParmCodes(2)=='f',obj.sigma,0.95*obj.sigma);
            obj.ResetParms([NewMean NewSigma]);
        end
        
        function []=ReInit(obj)
            ReInit@RNGamma(obj);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
    end  % methods
    
end  % class RNGammaMS




