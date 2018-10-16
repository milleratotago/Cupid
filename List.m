classdef List < dDiscrete
    % List(Xs,Ps) is the distribution of any discrete set of X values (which should be sorted in increasing order).
    % Ps is a vector of the probabilities associated with each X, and these should sum to 1.0.
    
    methods
        
        function obj=List(varargin)
            obj=obj@dDiscrete('List');
            obj.ParmTypes = '';
            obj.DefaultParmCodes = '';
            obj.NDistParms = 0;
            
            obj.PDFNearlyZero = 2*eps(0);
            
            switch nargin
                case 0
                case 2
                    ResetParms(obj,varargin{:});
                otherwise
                    ME = MException('List:Constructor', ...
                        'List constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,Xs,Ps)
            ClearBeforeResetParmsD(obj);
            obj.DiscreteX = Xs;
            obj.DiscretePDF = Ps;
            ReInit(obj);
        end
        
%        function []=PerturbParms(obj,ParmCodes)
%            ReInit(obj);
%        end
        
        function []=ReInit(obj)
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.Initialized = true;
            TrimTables(obj,eps(0),1);
            SetBinEdges(obj);
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~) %#ok<INUSL>
%             Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.Bounded2Real(0,1,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~) %#ok<INUSL>
%             Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2Bounded(0,1,Reals(2))];
        end
               
    end  % methods
    
end  % class List


