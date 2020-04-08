classdef RecinormalInv < Recinormal
    % Distribution of Y = 1/Z, where Z is normally distributed.
    % The parameters of this version are 1/mu_z and sqrt(1/sigma_z),
    % which are more convenient units for parameter estimation in some cases.

    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        OneOverMu, OneOverSigma
    end
    
    methods
        
        function obj=RecinormalInv(OneOverMu,OneOverSigma)   % Constructor
            obj=obj@Recinormal;
            obj.FamilyName = 'RecinormalInv';
            obj.ParmNames{1} = 'OneOverMu';
            obj.ParmNames{2} = 'OneOverSigma';
            ResetParms(obj,[OneOverMu,OneOverSigma]);
        end

        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newmu    = ifelse(ParmCodes(1)=='f', obj.OneOverMu,    1.01*obj.OneOverMu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.OneOverSigma, 0.99*obj.OneOverSigma);
            obj.ResetParms([newmu newsigma]);
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.OneOverMu = newparmvalues(1);
            obj.OneOverSigma = newparmvalues(2);
            obj.mu = 1 / obj.OneOverMu;
            obj.sigma = 1 / obj.OneOverSigma;
            ReInit(obj);
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.OneOverMu) ',' num2str(obj.OneOverSigma) ')'];
        end
        
    end % methods

end % classdef
