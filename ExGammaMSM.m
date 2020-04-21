classdef ExGammaMSM < ExGamma
    % ExGammaMSM distribution (sum of gamma and exponential) with parameters gamma mean, gamma SD, exponential mean.
    % Constraint to avoid numerical errors: muG >= sigmaG
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        muG, sigmaG, muE
    end
    
    methods
        
        function obj=ExGammaMSM(varargin)
            obj=obj@ExGamma;
            obj.FamilyName = 'ExGammaMSM';
            obj.ParmNames{1} = 'muG';
            obj.ParmNames{2} = 'sigmaG';
            obj.ParmNames{3} = 'muE';
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExGammaMSM:Constructor', ...
                        'ExGammaMSM constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.muG = newparmvalues(1);
            obj.sigmaG = newparmvalues(2);
            obj.muE = newparmvalues(3);
            tempKG = (obj.muG/obj.sigmaG)^2;
            tempRateG = tempKG / obj.muG;
            tempRateE = 1 / obj.muE;
            HoldNameBuilding = obj.NameBuilding;
            obj.NameBuilding = false;
            ResetParms@ExGamma(obj,[tempKG tempRateG tempRateE]);
            obj.NameBuilding = HoldNameBuilding;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, just make everything a little larger.
            newMuG    = ifelse(ParmCodes(1)=='f', obj.muG,   1.05*obj.muG);
            newSigmaG = ifelse(ParmCodes(2)=='f', obj.sigmaG,1.02*obj.sigmaG);
            newMuE  = ifelse(ParmCodes(3)=='f', obj.muE,     1.03*obj.muE);
            obj.ResetParms([newMuG newSigmaG newMuE]);
        end
        
        function []=ReInit(obj)
            ReInit@ExGamma(obj);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            muG = Parms(1);
            Reals = [NumTrans.GT2Real(eps,muG) NumTrans.Bounded2Real(eps,muG,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            muG = NumTrans.Real2GT(eps,Reals(1));
            Parms = [muG NumTrans.Real2Bounded(eps,muG,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
    end  % methods
    
end  % class ExGammaMSM
