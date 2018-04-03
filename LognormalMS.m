classdef LognormalMS < Lognormal
    % LognormalMS(postmu,postsigma):  parameters are the final mean & sd, after exp()
    
    properties(SetAccess = protected)
        postmu, postsigma
    end
    
    methods
        
        function obj=LognormalMS(varargin)   % Constructor
            obj=obj@Lognormal;
            obj.ThisFamilyName = 'LognormalMS';
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
            CheckBeforeResetParms(obj,newparmvalues);
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
            assert(obj.postmu>0,'LognormalMS postmu must be > 0.');
            assert(obj.postsigma>0,'LognormalMS postsigma must be > 0.');
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
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
    end  % methods
    
end  % class LognormalMS

