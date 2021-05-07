classdef Uquad < dContinuous
    % Uquad(min,max) quadratic distribution, with minimum at center.
    % See https://en.wikipedia.org/wiki/U-quadratic_distribution
    
    properties(SetAccess = protected)
        min, max
        alpha, beta  % derived from min, max
    end
    
    methods
        
        function obj=Uquad(varargin)
            obj=obj@dContinuous('Uquad');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Uquad:Constructor', ...
                        'Uquad constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.min = newparmvalues(1);
            obj.max = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, make the bounds a bit wider.
            OldLower = obj.LowerBound;
            OldUpper = obj.UpperBound;
            BoundShift = 0.051 * (OldUpper - OldLower);
            NewLower = ifelse(ParmCodes(1)=='f',obj.LowerBound,OldLower-BoundShift);
            NewUpper = ifelse(ParmCodes(2)=='f',obj.UpperBound,OldUpper+BoundShift);
            obj.ResetParms([NewLower NewUpper]);
        end
        
        function []=ReInit(obj)
            assert(obj.min<obj.max,'Uquad parameters must be min<max.');
            obj.beta = (obj.min + obj.max) / 2;
            obj.alpha = 12 / (obj.max-obj.min)^3;
            obj.LowerBound = obj.min;
            obj.UpperBound = obj.max;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % It is tempting to enforce the restriction that Parms(2) has to be greater than Parms(1), like this:
            % Reals = [Parms(1) NumTrans.GT2Real(Parms(1),Parms(2))];
            % Unfortunately, this creates problems if you want to estimate Parms(1) while holding Parms(2) fixed.
            % It may be better to use a different parameterization, e.g. Parms(1) = beta & Parms(2) = width,
            % so that these parameters are independent, but this has not been implemented.
            Reals = Parms;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = Reals;
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            InBounds = (X >= obj.LowerBound) & (X <= obj.UpperBound);
            thispdf(InBounds) = obj.alpha * (X(InBounds) - obj.beta).^2;
        end% 
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            thiscdf(X>obj.UpperBound) = 1;
            InBounds = (X >= obj.LowerBound) & (X <= obj.UpperBound);
            thiscdf(InBounds) = obj.alpha/3 * ( (X(InBounds) - obj.beta).^3 + (obj.beta-obj.min)^3 );
        end

%        function thisval=InverseCDF(obj,P)  started with symbolic at the end
%            assert(obj.Initialized,UninitializedError(obj));
%            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
%            thisval = zeros(size(P));
%        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.beta;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 3 * (obj.max-obj.min)^2 / 20;
        end
        
        function thisval=RawSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0;
        end
        
%         function thisval=Kurtosis(obj)
%             assert(obj.Initialized,UninitializedError(obj));
%             thisval = 3 / 112 * (obj.max - obj.min)^4;
%         end
        
%        function thisval=Random(obj,varargin)  NEWJEFF
%            assert(obj.Initialized,UninitializedError(obj));
%        end
        
    end  % methods
    
end  % class Uquad

% Symbolic computations for determining Uquad inversecdf
% syms p x alpha beta
% CDFfn = p == alpha/3 * ( (x-beta)^3 + (beta-alpha)^3 );
% solution = solve(CDFfn,x)
