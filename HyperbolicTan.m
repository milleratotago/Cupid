classdef HyperbolicTan < dContinuous
    % Hyperbolic Tangent distribution with parameter Scale>0.
    %  See Strasburger (2002, Perception & Psychophysics.}
    
    properties(SetAccess = protected)
        Scale
    end
    
    methods
        
        function obj=HyperbolicTan(varargin)
            obj=obj@dContinuous('HyperbolicTan');
            obj.NDistParms = 1;
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('HyperbolicTan:Constructor', ...
                        'HyperbolicTan constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.Scale = newparmvalues(1);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newScale = ifelse(ParmCodes(1)=='f', obj.Scale, 0.9*obj.Scale);
            obj.ResetParms(newScale);
        end
        
        function []=ReInit(obj)
            assert(obj.Scale>0,'HyperbolicTan Scale parameter must be > 0.');
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.GT2Real(eps,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(eps,Reals(1));
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            ExpPlus = exp(obj.Scale*X(InBounds));
            ExpMinus = 1 ./ ExpPlus;
            thispdf(InBounds) = 4 * obj.Scale ./ (ExpPlus + ExpMinus).^2;
            % for i=1:numel(X)
            %     if (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
            %         ExpPlus = exp(obj.Scale*X(i));
            %         ExpMinus = 1/ExpPlus;
            %         thispdf(i) = 4 * obj.Scale / (ExpPlus + ExpMinus)^2;
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(X>=obj.UpperBound) = 1;
            ExpPlus = exp(obj.Scale*X(InBounds));
            ExpMinus = 1 ./ ExpPlus;
            thiscdf(InBounds) = (ExpPlus-ExpMinus) ./ (ExpPlus+ExpMinus);
            % for i=1:numel(X)
            %     if X(i)<=obj.LowerBound
            %     elseif X(i)>=obj.UpperBound
            %         thiscdf(i) = 1;
            %     else
            %         ExpPlus = exp(obj.Scale*X(i));
            %         ExpMinus = 1/ExpPlus;
            %         thiscdf(i) = (ExpPlus - ExpMinus) / (ExpPlus+ExpMinus);
            %     end
            % end
        end
        
        function thisval=InverseCDF(obj,P)  % NWJEFF: Vectorize InverseCDFs. InBounds?
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            OMP = 1 - P;
            S = sqrt(OMP.*(1+P)) ./ OMP;
            thisval = log(S) / obj.Scale;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = log(2) / obj.Scale;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            EXSqr = (pi/obj.Scale)^2 / 12.0;
            thisval = EXSqr - Mean(obj)^2;
        end
        
        function thisval=Median(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0.5493 / obj.Scale;
        end
        
    end  % methods
    
end  % class HyperbolicTan

