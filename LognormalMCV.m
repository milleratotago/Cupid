classdef LognormalMCV < Lognormal
    % LognormalMCV(postmu,postsigmapct):  parameters are the final mean & CV=sd/mu, after exp()
    
    properties(SetAccess = protected)
        postmu, postcv
        postsigma
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=LognormalMCV(varargin)   % Constructor
            obj=obj@Lognormal;
            obj.FamilyName = 'LognormalMCV';
            obj.ParmNames{1} = 'postmu';
            obj.ParmNames{2} = 'postcv';
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('LognormalMCV:Constructor', ...
                        'LognormalMCV constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.postmu = newparmvalues(1);
            obj.postcv = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            NewMean = ifelse(ParmCodes(1)=='f',obj.postmu,1.05*obj.postmu);
            Newcv = ifelse(ParmCodes(2)=='f',obj.postsigma,0.95*obj.postcv);
            obj.ResetParms([NewMean Newcv]);
        end
        
        function []=ReInit(obj)
            assert(obj.postmu>0,'LognormalMCV postmu must be > 0.');
            assert(obj.postcv>0,'LognormalMCV postcv must be > 0.');
            obj.postsigma = obj.postcv * obj.postmu;
            obj.sigma = sqrt( log( obj.postsigma^2 / obj.postmu^2 + 1) );
            obj.mu = log( obj.postmu/exp(0.5*obj.sigma^2) );
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function s = EstMom(obj,TargetVals,varargin)
            if numel(varargin)<1
                ParmCodes = obj.DefaultParmCodes;
            else
                ParmCodes = varargin{1};
            end
            if (numel(TargetVals)==1)&&(ParmCodes(1)=='r')
                obj.ResetParms([TargetVals(1), obj.postcv]);
                s = obj.StringName;
            elseif (numel(TargetVals)==2)&&(ParmCodes(1)=='r')&&(ParmCodes(2)=='r')
                WantCV = sqrt(TargetVals(2))/TargetVals(1);
                obj.ResetParms([TargetVals(1), WantCV]);
                s = obj.StringName;
            else
                s = EstMom@dGeneric(obj,TargetVals,varargin{:});
            end
        end

    end  % methods
    
end  % class LognormalMCV

