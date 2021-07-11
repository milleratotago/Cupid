classdef Chi < dContinuous
    % Chi distribution with integer degrees of freedom df > 0.
    
    properties(SetAccess = protected)
        df,  % distribution parameter degrees of freedom
        dfm1, Halfdf, Halfdfm1, GammaHalfdf, LnGammaHalfdf, Inv2toHalfdfm1GammaHalfdf
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            % Avoid mapping with integer parms
            Reals = Parms;  % [NumTrans.GT2Real(1,Parms(1))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = Reals; % [NumTrans.Real2GT(1,Reals(1))];
        end
        
    end
    
    methods
        
        function obj=Chi(varargin)   % Constructor
            obj=obj@dContinuous('Chi');
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            obj.NDistParms = 1;
            % Handle constructor calls with different numbers of parameters.
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('Chi:Constructor', ...
                        'Too many arguments passed to Chi constructor.');
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
            obj.dfm1 = obj.df - 1;
            obj.Halfdf = obj.df / 2 ;
            obj.Halfdfm1 = obj.Halfdf - 1;
            obj.GammaHalfdf = gamma(obj.Halfdf);
            obj.LnGammaHalfdf = log(obj.GammaHalfdf);
            obj.Inv2toHalfdfm1GammaHalfdf = 1 / ( 2^obj.Halfdfm1 * obj.GammaHalfdf );
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
            % JKB Vol 1 Eqn 18.5
            thispdf(InBounds) = obj.Inv2toHalfdfm1GammaHalfdf * exp(-X(InBounds).^2/2) .* X(InBounds).^obj.dfm1;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            % JKB Vol 1 Eqn 18.3
            thiscdf(InBounds) = gammainc(X(InBounds).^2/2,obj.Halfdf);  % Was GammaIntegral(X^2/2,obj.Halfdf,obj.LnGammaHalfdf)
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % JKB Vol 1 Eqn 18.14
            thisval = sqrt(2) * gamma((obj.df+1)/2) / obj.GammaHalfdf;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % JKB Vol 1 Eqn 18.14
            thisval = obj.df - 2 * ( gamma((obj.df+1)/2) / obj.GammaHalfdf )^2;
        end
        
        function thisval=RawMoment(obj,I)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % JKB Vol 1 Eqn 18.13
            thisval = 2^(I/2) * gamma((obj.df+I)/2) / obj.GammaHalfdf;
        end
        
    end  % methods
    
end  % class Chi

