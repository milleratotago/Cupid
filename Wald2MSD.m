classdef Wald2MSD < Wald  % NEWJEFF: Not documented or unit-tested
    % 2-parameter version of Wald distribution (mu,sigma,barrierA) with sigma=1
    % The two parameters of this version of the distribution are its mean and SD.
    
    properties(SetAccess = protected)
        mean, sd
    end

    properties(SetAccess = public)
        minmean, minsd
    end

    methods
        
        function obj=Wald2MSD(varargin)
            obj=obj@Wald;
            obj.FamilyName = 'Wald2MSD';
            obj.ParmNames{1} = 'mean';
            obj.ParmNames{2} = 'sd';
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.minmean = eps;
            obj.minsd = eps;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Wald2MSD:Constructor', ...
                        'Wald2MSD constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mean = newparmvalues(1);
            obj.sd = newparmvalues(2);
            assert(obj.mean>0,'Wald2MSD mean must be > 0.');
            assert(obj.sd>0,'Wald2MSD sd must be > 0.');
            passMu = sqrt(obj.mean/obj.sd^2);
            passbarrierA = passMu * obj.mean;
            obj.NDistParms = 3;  % Pretend there are 3 parameters for inherited Wald
            ResetParms@Wald(obj,[passMu 1 passbarrierA]);
            obj.NDistParms = 2;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newmean    = ifelse(ParmCodes(1)=='f', obj.mean,      1.05*obj.mean);
            newsd     = ifelse(ParmCodes(2)=='f', obj.sd,1.05*obj.sd);
            obj.ResetParms([newmean newsd]);
        end
        
        function []=BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.mean) ',' num2str(obj.sd) ')'];
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(obj.minmean,Parms(1)) NumTrans.GT2Real(obj.minsd,Parms(2))] ;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(obj.minmean,Reals(1)) NumTrans.Real2GT(obj.minsd,Reals(2))];
        end
        
    end  % methods
    
end  % class Wald2MSD


