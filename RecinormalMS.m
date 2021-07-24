classdef RecinormalMS < Recinormal
    % Distribution of Y = 1/Z, where Z is normally distributed.
    % This version is parameterized in terms of finalmu, finalsigma,
    % which are the mean and standard deviation of Y.
    % Note that parameter estimation with this version is extremely slow.
    % This is because for each new (finalmu, finalsigma) pair that is considered,
    % a numerical search must be carried out to find the appropriate mu,sigma for Z.

    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        finalmu, finalsigma
        utilityRecinormal
    end
    
    methods
        
        function obj=RecinormalMS(finalmu,finalsigma)   % Constructor
            obj=obj@Recinormal;
            obj.utilityRecinormal = Recinormal(.003,.0002);
            obj.FamilyName = 'RecinormalMS';
            ResetParms(obj,[finalmu, finalsigma]);
        end

        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newmu    = ifelse(ParmCodes(1)=='f', obj.finalmu,    1.01*obj.finalmu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.finalsigma, 0.99*obj.finalsigma);
            obj.ResetParms([newmu newsigma]);
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.finalmu = newparmvalues(1);
            obj.finalsigma = newparmvalues(2);
            newparmvalues(2) = newparmvalues(2)^2;
            % Infinite recursion problem if I use obj.EstMom directly.
            obj.utilityRecinormal.EstMom(newparmvalues);
            obj.mu = obj.utilityRecinormal.mu;
            obj.sigma = obj.utilityRecinormal.sigma;
            ReInit(obj);
        end

        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.finalmu) ',' num2str(obj.finalsigma) ')'];
        end
        
        function parms = StartParmsMLE(obj,X)
            est_mu = mean(X);
            est_sigma = std(X);
            parms = [est_mu, est_sigma];
        end
        
    end % methods

end % classdef
