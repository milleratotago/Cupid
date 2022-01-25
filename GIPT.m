classdef GIPT < dContinuous % dEither  % Discrete not implemented, like AttainP
    % GIPT(X,P): P is a distribution of 0-1 values.  X is any (continuous) distribution.
    % GIPT is the distribution formed by transforming the P values via X.InverseCDF(P).
    
    properties(SetAccess = protected)
        X, P
    end
    
    properties(SetAccess = private)
        NXParms   % Holds the number of parameters in the X distribution.
                  % These appear at the start of parameter lists.
    end
    
    methods
        
        function obj=GIPT(varargin)
            obj=obj@dContinuous('GIPT');
            switch nargin
                case 0
                case 2
                    Setup(obj,varargin(:));
                    ResetParms(obj,[obj.X.ParmValues obj.P.ParmValues]);
                otherwise
                    ME = MException('GIPT:Constructor', ...
                        'GIPT constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function Setup(obj,s)
            obj.X = s{1};
            obj.P = s{2};
            obj.DistType = obj.X.DistType;
        end
        
        function BuildMyName(obj)
            obj.StringName = ['GIPT(' obj.X.StringName ',' obj.P.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Initialized = false;
            obj.X.ResetParms(newparmvalues(1:obj.X.NDistParms));
            obj.P.ResetParms(newparmvalues(obj.X.NDistParms+1:end));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.X.PerturbParms(ParmCodes(1:obj.X.NDistParms));
            obj.P.PerturbParms(ParmCodes(obj.X.NDistParms+1:end));
            obj.ResetParms([obj.X.ParmValues obj.P.ParmValues]);
        end
        
        function []=ReInit(obj)
            obj.NXParms = obj.X.NDistParms;
            obj.NDistParms = obj.P.NDistParms + obj.X.NDistParms;
            % By default, parameters of the P distribution are fixed
            % and those of the X distribution are allowed to vary.
            s = obj.P.DefaultParmCodes;
            s(1:end) = 'f';
            obj.DefaultParmCodes = [obj.X.DefaultParmCodes s];
            obj.Initialized = true;
            obj.DistType = obj.X.DistType;
            obj.NValues = obj.X.NValues;
            obj.LowerBound = obj.X.InverseCDF(obj.P.LowerBound);
            obj.UpperBound = obj.X.InverseCDF(obj.P.UpperBound);
            if obj.NameBuilding
                BuildMyName(obj);
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            p = obj.X.CDF(X(InBounds));
            thiscdf(InBounds) = obj.P.CDF(p);
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.X.ParmsToReals(Parms(1:obj.X.NDistParms)) obj.P.ParmsToReals(Parms(obj.X.NDistParms+1:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.X.RealsToParms(Reals(1:obj.X.NDistParms)) obj.P.RealsToParms(Reals(obj.X.NDistParms+1:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.X.ParmValues obj.P.ParmValues];
        end
        
        % function thisval = nIthValue(iVal)
        %     thisval = 1;  % Dummy
        % end
        
        % function thisval = NearestLegal(X)
        %     thisval = 1;  % Dummy
        % end
        
        % function thisval = LegalValue(X)
        %     thisval = false;  % Dummy
        % end
        
        function x = PtoX(obj,p)
            % Convert p values from the P distribution into x values
            % relative to the X distribution.
            x = obj.X.InverseCDF(p);
        end
        
        function x = XsToPlot(obj)
            x = obj.PtoX(obj.P.XsToPlot);
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            p = obj.P.Random(varargin{:});
            thisval = obj.PtoX(p);
        end
        
    end  % methods
    
end  % class GIPT


