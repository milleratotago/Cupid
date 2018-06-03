classdef TriangularCW < Triangular
    % TriangularCW(center,width) distribution, with peak at center.
    
    properties(SetAccess = protected)
        width
        minwidth
    end
    
    methods
        
        function obj=TriangularCW(varargin)
            obj=obj@Triangular;
            obj.FamilyName = 'TriangularCW';
            obj.ParmNames{1} = 'center';
            obj.ParmNames{2} = 'width';
            obj.minwidth = sqrt(eps);
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('TriangularCW:Constructor', ...
                        'TriangularCW constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
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
            assert(obj.width>=obj.minwidth,'TriangularCW width must be > 0.');
            obj.min = obj.center - obj.width/2;
            obj.max = obj.center + obj.width/2;
            ReInit@Triangular(obj);
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(obj.minwidth,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(obj.minwidth,Reals(2))];
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            TriangularCW1 = rand(varargin{:});
            TriangularCW2 = rand(varargin{:});
            thisval = (TriangularCW1 - TriangularCW2) * (obj.center - obj.min) + obj.center;
        end
        
    end  % methods
    
end  % class TriangularCW


