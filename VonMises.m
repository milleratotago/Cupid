classdef VonMises < dContinuous
    % VonMises(Loc,Concentration>0)
    % Support is on interval Loc-pi to Loc+pi
    % https://en.wikipedia.org/wiki/Von_Mises_distribution
    
    properties(SetAccess = protected)
        Loc, Concentration
        BesselEps, Bessel0Scale, TwoPi, TwoPiBessel0Scale
    end
    
    methods
        
        function obj=VonMises(varargin)
            obj=obj@dContinuous('VonMises');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.BesselEps = 1e-10;
            obj.TwoPi = 2*pi;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('VonMises:Constructor', ...
                        'VonMises constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.Loc = newparmvalues(1);
            obj.Concentration = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newloc   = ifelse(ParmCodes(1)=='f', obj.Loc, 0.9*obj.Loc);
            newscale = ifelse(ParmCodes(2)=='f', obj.Concentration,1.1*obj.Concentration);
            obj.ResetParms([newloc newscale]);
        end
        
        function []=ReInit(obj)
            assert(obj.Concentration>0,'VonMises Concentration must be > 0.');
            obj.Bessel0Scale = besseli(0,obj.Concentration);  % pascal was ModBessI0(obj.Concentration);
            obj.TwoPiBessel0Scale = obj.TwoPi * obj.Bessel0Scale;
            obj.LowerBound = obj.Loc - pi;
            obj.UpperBound = obj.Loc + pi;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(0,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(0,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = exp(obj.Concentration*cos(X(InBounds)-obj.Loc)) / obj.TwoPiBessel0Scale;
            % for i=1:numel(X)
            %     if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %         thispdf(i) = exp(obj.Concentration*cos(X(i)-obj.Loc)) / obj.TwoPiBessel0Scale;
            %     end
            % end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.Loc;
        end
        
    end  % methods
    
end  % class VonMises


