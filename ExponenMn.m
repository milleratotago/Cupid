classdef ExponenMn < Exponential
    % ExponenMn distribution: exponential with parameter = mean, not rate.
    
    properties(SetAccess = protected)
        expmean
    end
    
    methods
        
        function obj=ExponenMn(varargin)
            obj=obj@Exponential;
            obj.ThisFamilyName = 'ExponenMn';
            obj.ParmNames{1} = 'expmean';
            if nargin==1
                ResetParms(obj,varargin{1});
            elseif nargin>1
                ME = MException('Exponential:Constructor', ...
                    'Too many arguments passed to Exponential constructor.');
                throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.expmean = newparmvalues(1);
            obj.rate = 1/obj.expmean;
            assert(obj.expmean>0,'ExponenMn mean must be > 0.');
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newmean = ifelse(ParmCodes(1)=='f', obj.expmean, 0.9*obj.expmean);
            obj.ResetParms(newmean);
        end
        
        %        function parmvals = ParmValues(obj,varargin)
        %            parmvals = obj.expmean;
        %        end
        
        function s=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            meanObs = mean(Observations);
            ResetParms(obj,meanObs);
            BuildMyName(obj);
            s=obj.StringName;
        end
        
    end
    
end

