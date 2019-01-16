classdef InverseGaussian < Wald  % NEWJEFF: No documentation or unit tests.
    % 2-parameter version of Wald distribution with parameters IGmu and lambda

    properties
        IGmu, % mu refers to the Wald mu, which is a different parameter
        lambda
    end
    
    methods
        
        function obj=InverseGaussian(varargin)
            obj=obj@Wald;
            obj.FamilyName = 'InverseGaussian';
            obj.ParmNames{1} = 'mu';
            obj.ParmNames{2} = 'lambda';
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('InverseGaussian:Constructor', ...
                        'InverseGaussian constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.IGmu = newparmvalues(1);
            obj.lambda = newparmvalues(2);
            assert(obj.IGmu>0,'InverseGaussian mu must be > 0.');
            assert(obj.lambda>0,'InverseGaussian lambda must be > 0.');
            obj.NDistParms = 3;  % Pretend there are 3 parameters for inherited Wald
            passMu = sqrt(obj.lambda)/obj.IGmu;
            passbarrierA = sqrt(obj.lambda);
            ResetParms@Wald(obj,[passMu 1 passbarrierA]);
            obj.NDistParms = 2;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newmu    = ifelse(ParmCodes(1)=='f', obj.IGmu,  1.05*obj.IGmu);
            newlam   = ifelse(ParmCodes(2)=='f', obj.lambda,1.05*obj.lambda);
            obj.ResetParms([newmu newlam]);
        end
        
        function []=BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.IGmu) ',' num2str(obj.lambda) ')'];
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(0,Parms(1)) NumTrans.GT2Real(0,Parms(2))] ;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(0,Reals(1)) NumTrans.Real2GT(0,Reals(2))];
        end
        
    end  % methods
    
end  % class InverseGaussian


