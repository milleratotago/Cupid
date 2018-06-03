classdef ChiSqNoncentral < dContinuous
    % Noncentral Chi^2 distribution with degrees of freedom df > 0 (must be an integer).
    % noncen (noncentrality) is the sum across the df chi-squared RVs of the squared
    % mean of the normal going into that chi-square (divided by its sd)
    
    properties(SetAccess = protected)
        df, noncen,
        Halfdf, ExpNegHalfNoncen, MeanOfEachNormal
    end
    
    methods
        
        function obj=ChiSqNoncentral(varargin)
            obj=obj@dContinuous('ChiSqNoncentral');
            obj.ParmTypes = 'fr';
            obj.DefaultParmCodes = 'fr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ChiSqNoncentral:Constructor', ...
                        'ChiSqNoncentral:Constructor requires 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.df = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.noncen = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newdf     = ifelse(ParmCodes(1)=='f', obj.df,      obj.df+1);
            newnoncen = ifelse(ParmCodes(2)=='f', obj.noncen,   1.1*obj.noncen);
            obj.ResetParms([ newdf newnoncen ]);
        end
        
        function []=ReInit(obj)
            assert((obj.df>0)&&(iswholenumber(obj.df)),'ChiSqNoncentral df must be an integer > 0.');
            assert(obj.noncen>=0,'ChiSqNoncentral noncen must be >= 0.');
            obj.Halfdf = obj.df / 2 ;
            obj.ExpNegHalfNoncen = exp(-obj.noncen/2);
            obj.MeanOfEachNormal = sqrt(obj.noncen/obj.df);
            obj.LowerBound = obj.XNearlyZero;
            obj.UpperBound = 4000;
            obj.Initialized = true;  % Needed so that CDF can be called
            %while CDF(obj,2*obj.LowerBound) < obj.CDFNearlyZero
            %    obj.LowerBound = obj.LowerBound * 2;
            %end
            %while CDF(obj,obj.UpperBound) < obj.CDFNearlyOne
            %    obj.UpperBound = obj.UpperBound * 2;
            %end
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.GT2Real(0,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2GT(0,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            % Wikipedia:
            thispdf(InBounds) = 0.5*exp(-(X(InBounds)+obj.noncen)/2).*(X(InBounds)/obj.noncen).^(obj.df/4-0.5).*besseli(obj.Halfdf-1,sqrt(obj.noncen*X(InBounds)));
            % for iel=1:numel(X)
            %     if InBounds(iel)
            %          % Speed up here by computing & saving coeffs in init.
            %          thispdf(iel) = 0.5*exp(-(X(iel)+obj.noncen)/2)*(X(iel)/obj.noncen)^(obj.df/4-0.5)*besseli(obj.Halfdf-1,sqrt(obj.noncen*X(iel)));
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)   % JKB Vol 2 Eqn 29.2
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel = 1:numel(X)
                if InBounds(iel)
                    sum = 0;
                    j = -1;
                    jfac = 1;
                    term = 99;
                    % It should be possible to speed this up by pre-computing
                    % & saving coeff = numerator / denom  in init (or perhaps only when
                    % CDF is called).
                    while abs(term) > 1e-10
                        j = j + 1;
                        if j > 0
                            jfac = jfac*j;
                            nextgamma = nextgamma*(obj.Halfdf+j-1);
                        else
                            nextgamma = gamma(obj.Halfdf);
                        end
                        numerator = (0.5*obj.noncen)^j / jfac;
                        denom = 2^(obj.Halfdf+j) * nextgamma;
                        coeff = numerator / denom;
                        thisintegral = integral(@(y) y.^(obj.Halfdf+j-1) .* exp(-y/2), 0, X(iel));
                        term = coeff * thisintegral;
                        sum = sum + term;
                    end
                    thiscdf(iel) = obj.ExpNegHalfNoncen * sum;
                end
            end
        end
        
        function thisval=Mean(obj)
            % Wikipedia:
            thisval=obj.df + obj.noncen;
        end
        
        function thisval=Variance(obj)
            % Wikipedia:
            thisval=2*(obj.df+2*obj.noncen);
        end
        
        function thisval=RelSkewness(obj)
            % Wikipedia:
            thisval=2^(1.5)*(obj.df+3*obj.noncen)/(obj.df+2*obj.noncen)^1.5;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            try
                % Sum of df copies of random normal(MeanOfEachNormal,1) rvs:
                randnors = randn(varargin{:},obj.df) + obj.MeanOfEachNormal;
                randnors = randnors.^2;
                % sum only produces output with the desired number of dimensions if summing is across the last dimension.
                thisval = sum(randnors,ndims(randnors));
            catch
                % ... when obj.df is not an integer
                thisval = Random@dGeneric(obj,varargin{:});
            end
        end
        
    end  % methods
    
end  % class ChiSqNoncentral

