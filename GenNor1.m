classdef GenNor1 < dContinuous
    % This is version 1 of https://en.wikipedia.org/wiki/Generalized_normal_distribution
    % "General Error Distribution",
    
    % Notes:
    % This family includes the normal distribution when \beta =2 (with mean \mu and variance \alpha^{2}/2 
    % and it includes the Laplace distribution when \beta =1.
    % As \beta \rightarrow \infty, the density converges pointwise to a uniform density on
    %  (\mu -\alpha ,\mu +\alpha ).
    
    % This code requires the gamma and incomplete gamma functions

    % The distribution may also be called the "Error Distribution" by Evans, Hasting, & Peacock (1993), p. 57.
    %  and "Subbotin's distribution" by Johnson, Kotz, & Balakrishnan, 1995, Vol 2, p. 195
    % GenNor1(mu,scale,1) = Laplace(mu,scale)
    % GenNor1(mu,scale,2) = Normal(mu,scale)
    
    properties(SetAccess = protected)
        Mu     % Location, Real
        Alpha  % scale, positive real
        Beta   % shape, positive real
        PDFmul % Constant multiplier used in PDF
    end
    
    methods
        
        function obj=GenNor1(varargin)
            obj=obj@dContinuous('GenNor1');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.SearchOptions.MaxFunEvals = 20000;
            obj.SearchOptions.MaxIter = 10000;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('GenNor1:Constructor', ...
                        'GenNor1 constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.Mu = newparmvalues(1);
            obj.Alpha = newparmvalues(2);
            obj.Beta = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            NewMu  = ifelse(ParmCodes(1)=='f', obj.Mu,obj.Mu + 0.5);
            NewAlpha = ifelse(ParmCodes(2)=='f', obj.Alpha,1.1*obj.Alpha);
            NewBeta = ifelse(ParmCodes(3)=='f', obj.Beta,1.1*obj.Beta);
            obj.ResetParms([NewMu NewAlpha NewBeta]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = false;
            assert(obj.Alpha>0,'GenNor1 Alpha must be > 0.');
            assert(obj.Beta>0,'GenNor1 Beta must be > 0.');
            obj.PDFmul = obj.Beta / (2*obj.Alpha*gamma(1/obj.Beta)); 

            obj.LowerBound = obj.Mu - 100*obj.Alpha;
            obj.UpperBound = obj.Mu + 100*obj.Alpha;
            obj.Initialized = true;
            obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
            obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            Z = (abs(X - obj.Mu)/obj.Alpha).^obj.Beta;
            thispdf = obj.PDFmul * exp(-Z);
            thispdf(X<obj.LowerBound) = 0;
            thispdf(X>obj.UpperBound) = 0;
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            sgn = sign(X-obj.Mu);
            % Note: The following line looks a little different than what is shown in Wikipedia
            % because it and MATLAB define gammainc differently by a factor of gamma(1/obj.Beta).
            frac = gammainc((abs(X-obj.Mu)/obj.Alpha).^obj.Beta, 1/obj.Beta)/2;
            thiscdf = 0.5 + sgn.*frac;
            thiscdf(X<obj.LowerBound) = 0;
            thiscdf(X>obj.UpperBound) = 1;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.Mu;
        end
        
        function thisval=Median(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.Mu;
        end
        
        function thisval=Variance(obj)
            thisval = obj.Alpha^2*gamma(3/obj.Beta) / gamma(1/obj.Beta);
        end
        
        function thisval=RawSkewness(obj)
            thisval = 0;
        end
        
        function thisval=Kurtosis(obj)
            thisval = gamma(5/obj.Beta)*gamma(1/obj.Beta) / gamma(3/obj.Beta)^2;
        end
        
    end  % methods
    
    
end  % class GenNor1
