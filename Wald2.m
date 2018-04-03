classdef Wald2 < Wald
    % 2-parameter version of Wald distribution (mu,sigma,barrierA) with sigma=1
    
    methods
        
        function obj=Wald2(varargin)
            obj=obj@Wald;
            obj.ThisFamilyName = 'Wald2';
            obj.ParmNames{2} = 'barrierA';
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Wald2:Constructor', ...
                        'Wald2 constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            passMu = newparmvalues(1);
            passbarrierA = newparmvalues(2);
            assert(passMu>0,'Wald2 mu must be > 0.');
            assert(passbarrierA>0,'Wald2 barrierA must be > 0.');
            obj.NDistParms = 3;  % Pretend there are 3 parameters for inherited Wald
            ResetParms@Wald(obj,[passMu 1 passbarrierA]);
            obj.NDistParms = 2;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newmu    = ifelse(ParmCodes(1)=='f', obj.mu,      1.05*obj.mu);
            newA     = ifelse(ParmCodes(2)=='f', obj.barrierA,1.05*obj.barrierA);
            obj.ResetParms([newmu newA]);
        end
        
        function []=BuildMyName(obj)
            obj.StringName = [obj.ThisFamilyName '(' num2str(obj.mu) ',' num2str(obj.barrierA) ')'];
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(0,Parms(1)) NumTrans.GT2Real(0,Parms(2))] ;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(0,Reals(1)) NumTrans.Real2GT(0,Reals(2))];
        end
        
    end  % methods
    
end  % class Wald2


