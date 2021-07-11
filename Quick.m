classdef Quick < dContinuous
    
    % Quick(Scale>0,Shape>0) distribution defined by Quick (1974, Kybernetics).
    %  My reference was Strasburger (2002, Perception & Psychophysics).
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        Scale, Shape, LnTwo
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(0,Parms(1)) NumTrans.GT2Real(0,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(0,Reals(1)) NumTrans.Real2GT(0,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=Quick(varargin)   % Constructor
            obj=obj@dContinuous('Quick');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.LnTwo = log(2);
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Quick:Constructor', ...
                        'Quick constructor requires 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newScale = ifelse(ParmCodes(1)=='f', obj.Scale, 1.1*obj.Scale);
            newShape = ifelse(ParmCodes(2)=='f', obj.Shape, 1.1*obj.Shape);
            obj.ResetParms([newScale newShape]);
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Scale = newparmvalues(1);
            obj.Shape = newparmvalues(2);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            assert(obj.Scale>0,'Quick Scale parameter must be > 0.');
            assert(obj.Shape>0,'Quick Shape parameter must be > 0.');
            obj.LowerBound = 0;
            obj.UpperBound =  1;
            obj.Initialized = true;
            while CDF(obj,obj.UpperBound) < obj.CDFNearlyOne
                obj.UpperBound = 2 * obj.UpperBound;
            end
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            y = X(InBounds) / obj.Scale;
            yp = y.^obj.Shape;
            thispdf(InBounds) =  2.^(-yp) .* yp * obj.Shape * obj.LnTwo ./ X(InBounds);
            %            for i=1:numel(X)
            %                thispdf(i) = 2^(-yp(i)) * yp(i) * obj.Shape * obj.LnTwo / X(i);
            %            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            y = X(InBounds) / obj.Scale;
            yp = y.^obj.Shape;
            TwoToYP = 2.^yp;
            thiscdf(InBounds) = 1 - 1 ./ TwoToYP;
            %            for i=1:numel(X)
            %                TwoToYP = 2^yp(i);
            %                thiscdf(i) = 1 - 1 / TwoToYP;
            %             end
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            for iel=1:numel(P)
                if InBounds(iel)
                    ExpStuff = log(1-P(iel)) / obj.LnTwo;
                    ExpStuff = log(-ExpStuff) / obj.Shape;
                    ExpStuff = exp(ExpStuff);
                    thisval(iel) = obj.Scale * ExpStuff;
                end
            end
        end
        
        function thisval=RawMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            P1 = obj.Scale^I;
            P2 = obj.LnTwo^(-I/obj.Shape);
            G  = gamma((I+obj.Shape)/obj.Shape);
            thisval = P1 * P2 * G;
        end
        
    end  % methods
    
end  % class Quick


