classdef ExGammaMSM < ExGamma
    % ExGammaMSM distribution (sum of gamma and exponential) with parameters gamma mean, gamma SD, exponential mean.
    % I once implemented the constraint muG >= sigmaG in an effort to avoid numerical errors,
    % but this is too restrictive so I have now relaxed it to muSigmaFac*muG >= sigmaG
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        muG, sigmaG, muE
    end
    
    properties(SetAccess = public)
        muSigmaFac
    end
    
    methods (Static)

        function parms = MomsToParms(GammaMean,GammaVar,ExpMean)
            % Return values of 3 distribution parameters yielding specified
            % mean and variance of normal and mean of exponential.
            % Used with ExHelpStartParmsMLE
            parms = zeros(3,1);
            parms(1) = GammaMean;
            parms(2) = sqrt(GammaVar);
            parms(3) = ExpMean;
        end

    end % methods (Static)

    methods
        
        function obj=ExGammaMSM(varargin)
            obj=obj@ExGamma;
            obj.FamilyName = 'ExGammaMSM';
            obj.ParmNames{1} = 'muG';
            obj.ParmNames{2} = 'sigmaG';
            obj.ParmNames{3} = 'muE';
            obj.muSigmaFac = 1.5;
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
            temp_g_shape = (obj.muG/obj.sigmaG)^2;
            temp_g_scale = obj.muG / temp_g_shape;
            HoldNameBuilding = obj.NameBuilding;
            obj.NameBuilding = false;
            ResetParms@ExGamma(obj,[temp_g_shape temp_g_scale obj.muE]);
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
        
%         function []=ReInit(obj)
%             ReInit@ExGamma(obj);
%             if (obj.NameBuilding)
%                 BuildMyName(obj);
%             end
%         end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % NEWJEFF: Here I assume muG>sigmaG but is that ever false?
            tmpmuG = Parms(1);
            sigmaGmax = obj.muSigmaFac * tmpmuG;
            if sigmaGmax < Parms(2)
                warning('obj.muSigmaFac*muG<sigmaG');
            end
            Reals = [NumTrans.GT2Real(eps,tmpmuG) NumTrans.Bounded2Real(eps,sigmaGmax,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            tmpmuG = NumTrans.Real2GT(eps,Reals(1));
            sigmaGmax = obj.muSigmaFac * tmpmuG;
            Parms = [tmpmuG NumTrans.Real2Bounded(eps,sigmaGmax,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
    end  % methods
    
end  % class ExGammaMSM
