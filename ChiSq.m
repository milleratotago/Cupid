classdef ChiSq < dContinuous
    % Chi^2 distribution with df degrees of freedom > 0.
    
    properties(SetAccess = protected)
        df,  % distribution parameter degrees of freedom
        Halfdf, Halfdfm1, GammaHalfdf, LnGammaHalfdf, Inv2toHalfdfGammaHalfdf
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = NumTrans.GT2Real(1,Parms(1));
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = NumTrans.Real2GT(1,Reals(1));
        end
        
    end
    
    methods
        
        function obj=ChiSq(varargin)   % Constructor
            obj=obj@dContinuous('ChiSq');
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            obj.NDistParms = 1;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('ChiSq:Constructor', ...
                        'Too many arguments passed to ChiSq constructor.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.df = VerifyIntegerGE(obj,1,newparmvalues(1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            if ~(ParmCodes(1) == 'f')
                obj.ResetParms(obj.df+1);
            end
        end
        
        function []=ReInit(obj)
            assert((obj.df>0)&&(iswholenumber(obj.df)),'ChiSq df must be an integer > 0.');
            obj.Halfdf = obj.df / 2 ;
            obj.Halfdfm1 = obj.Halfdf - 1;
            obj.GammaHalfdf = gamma(obj.Halfdf);
            obj.LnGammaHalfdf = log(obj.GammaHalfdf);
            obj.Inv2toHalfdfGammaHalfdf = 1 / ( 2^obj.Halfdf * obj.GammaHalfdf );
            obj.LowerBound = eps;
            obj.UpperBound = 4;
            obj.Initialized = true;  % Needed so that CDF can be called
            while CDF(obj,obj.UpperBound) < obj.CDFNearlyOne
                obj.UpperBound = obj.UpperBound * 2;
            end
            if obj.df > 10
            elseif obj.df > 1
                obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
                obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            else % df == 1 Begin
                obj.LowerBound = 1.0e-8;  % -9 gives incorrect PDFIntegral over lowest 10% of distribution
                obj.UpperBound = 18;
            end
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            % JKB Vol 1 Eqn 18.1
            thispdf(InBounds) = obj.Inv2toHalfdfGammaHalfdf * exp(-X(InBounds)/2) .* X(InBounds).^obj.Halfdfm1;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            % JKB Vol 1 Eqn 18.3
            thiscdf(InBounds) = gammainc(X(InBounds)/2,obj.Halfdf);  % Was GammaIntegral(X/2,obj.Halfdf,obj.LnGammaHalfdf)
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % JKB Vol 1 Eqn 18.10
            thisval = obj.df;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % JKB Vol 1 Eqn 18.10
            thisval = 2*obj.df;
        end
        
        function thisval=RawMoment(obj,I)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % JKB Vol 1 Eqn 18.8
            thisval = 1;
            for J = 1:I
                thisval = thisval * (obj.df + 2*(J-1));
            end
        end
        
        function thisval=MGF(obj,Theta) % From Johnson Kotz Balakrishnan
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            assert(abs(Theta)<0.5,'ChiSq MGF only defined for abs(Theta) < 0.5');
            thisval = (1-2*Theta)^(-obj.df/2);
        end
        
        function thisval=RelSkewness(obj) % From Mathematica
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 2 * sqrt(2) * sqrt(1/obj.df);
        end
        
        function thisval=Kurtosis(obj) % From Mathematica
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 3 + 12/ obj.df;
        end
        
    end  % methods
    
end  % class ChiSq

