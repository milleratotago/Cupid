classdef UniformCW < Uniform
    % UniformCW(center,width)
    
    properties(SetAccess = protected)
        center, width
        minwidth
    end
    
    methods
        
        function obj=UniformCW(varargin)
            obj=obj@Uniform;
            obj.ThisFamilyName = 'UniformCW';
            obj.ParmNames{1} = 'center';
            obj.ParmNames{2} = 'width';
            obj.minwidth = sqrt(eps);
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('UniformCW:Constructor', ...
                        'UniformCW:Constructor requires exactly 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.center = newparmvalues(1);
            obj.width = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            NewCenter = ifelse(ParmCodes(1)=='f',obj.center,obj.center+obj.width/10);
            NewWidth  = ifelse(ParmCodes(2)=='f',obj.width,1.1*obj.width);
            obj.ResetParms([NewCenter NewWidth]);
        end
        
        function []=ReInit(obj)
            assert(obj.width>0,'UniformCW must satisfy width>0.');
            obj.min = obj.center - obj.width/2;
            obj.max = obj.center + obj.width/2;
            ReInit@Uniform(obj);
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(obj.minwidth,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(obj.minwidth,Reals(2))];
        end
        
        function s=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            if numel(varargin)<1
                ParmCodes = obj.DefaultParmCodes;
            else
                ParmCodes = varargin{1};
            end
            minObs = min(Observations);
            maxObs = max(Observations);
            newcenter = ifelse(ParmCodes(1) == 'f',obj.center, (minObs+maxObs)/2);
            newwidth = ifelse(ParmCodes(2) == 'f',obj.width,maxObs-minObs);
            ResetParms(obj,[newcenter newwidth]);
            BuildMyName(obj);
            s=obj.StringName;
        end
        
    end  % methods
    
end  % class UniformCW

