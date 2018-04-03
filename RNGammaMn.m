classdef RNGammaMn < RNGamma
    % RNGammaMn(shape,mu) is an RNGamma specified by its shape parameter and its mean.

    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu
    end
    
    methods
        
        function obj=RNGammaMn(varargin)
            obj=obj@RNGamma;
            obj.ThisFamilyName = 'RNGammaMn';
            obj.ParmNames{1} = 'N';
            obj.ParmNames{2} = 'mu';
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('RNGammaMn:Constructor', ...
                        'RNGammaMn constructor must receive 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.N = newparmvalues(1);
            obj.mu = newparmvalues(2);
            assert(obj.N>0,'RNGammaMn N must be > 0.');
            assert(obj.mu>0,'RNGammaMn mu must be > 0.');
            tempRate = obj.N / obj.mu;
            HoldNameBuilding = obj.NameBuilding;
            obj.NameBuilding = false;
            ResetParms@RNGamma(obj,[obj.N tempRate]);
            obj.NameBuilding = HoldNameBuilding;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newN  = ifelse(ParmCodes(1)=='f',obj.N,1.1*obj.N);
            newMu = ifelse(ParmCodes(2)=='f',obj.mu,0.95*obj.mu);
            obj.ResetParms([newN newMu]);
        end
        
        function []=ReInit(obj)
            ReInit@RNGamma(obj);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
    end  % methods
    
end  % class RNGammaMn




