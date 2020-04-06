classdef LognormalMS < Lognormal
    % LognormalMS(postmu,postsigma):  parameters are the final mean & sd, after exp()
    
    properties(SetAccess = protected)
        postmu, postsigma
    end
    
    methods
        
        function obj=LognormalMS(varargin)   % Constructor
            obj=obj@Lognormal;
            obj.FamilyName = 'LognormalMS';
            obj.ParmNames{1} = 'postmu';
            obj.ParmNames{2} = 'postsigma';
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('LognormalMS:Constructor', ...
                        'LognormalMS constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.postmu = newparmvalues(1);
            obj.postsigma = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            NewMean = ifelse(ParmCodes(1)=='f',obj.postmu,1.05*obj.postmu);
            NewSigma = ifelse(ParmCodes(2)=='f',obj.postsigma,0.95*obj.postsigma);
            obj.ResetParms([NewMean NewSigma]);
        end
        
        function []=ReInit(obj)
            if obj.postmu<=0
                error('LognormalMS postmu must be > 0.');
            end
            if obj.postsigma<=0
                error('LognormalMS postsigma must be > 0.');
            end
            obj.sigma = sqrt( log( obj.postsigma^2 / obj.postmu^2 + 1) );
            obj.mu = log( obj.postmu/exp(0.5*obj.sigma^2) );
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(5*eps,Parms(1)) NumTrans.GT2Real(5*eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(5*eps,Reals(1)) NumTrans.Real2GT(5*eps,Reals(2))];
        end
        
        function s = EstMom(obj,TargetVals,varargin)
            s = EstMomMS(obj,TargetVals,varargin{:});
        end

    end  % methods
    
end  % class LognormalMS

