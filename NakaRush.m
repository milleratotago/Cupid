classdef NakaRush < dContinuous
    % Naka-Rushton(Scale>0):   distribution from psychophysics with CDF(x) = x^2/(1+x^2)
    
    properties(SetAccess = protected)
        Scale, ScaleSqr
    end
    
    methods
        
        function obj=NakaRush(varargin)   % Constructor
            obj=obj@dContinuous('NakaRush');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('NakaRush:Constructor', ...
                        'Too many arguments passed to NakaRush constructor.');
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
            if obj.Scale<5
                Factor = 1.1;
            else
                Factor = 0.9;
            end
            NewScale = ifelse(ParmCodes(1)=='f',obj.Scale,Factor*obj.Scale);
            obj.ResetParms(NewScale);
        end
        
        function []=ReInit(obj)
            assert(obj.Scale>0,'NakaRush Scale must be > 0.');
            obj.ScaleSqr = obj.Scale^2;
            obj.LowerBound = eps;
            obj.Initialized = true;
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);  % Boundary should match that in InverseCDF
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
            thispdf(InBounds) = 2 * X(InBounds) * obj.ScaleSqr ./ (1 + (obj.Scale*X(InBounds)).^2).^2 ;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            XSqr = (obj.Scale*X(InBounds)).^2;
            thiscdf(InBounds) = XSqr ./ (1+XSqr);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            Omp = (1 - P(InBounds)) * obj.ScaleSqr;
            thisval(InBounds) = sqrt( P(InBounds) ./ Omp );
            % for i=1:numel(thisval)
            %     if P(i) > obj.CDFNearlyOne
            %         thisval(i) = obj.UpperBound;
            %     else
            %         Omp = (1 - P(i)) * obj.ScaleSqr;
            %         thisval(i) = sqrt( P(i) / Omp );
            %     end
            % end
        end
        
        function thisval=Mean(obj)
            thisval = pi / (2 * obj.Scale);
        end
        
    end  % methods
    
end  % class NakaRush

