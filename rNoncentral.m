classdef rNoncentral < dContinuous
    % rNoncentral(SampleSize,TrueRho)
    
    % An older version using the Konishi (1978, Biometrika) approximation, etc, was saved in Cupid Storage.
    
    properties(SetAccess = protected)
        SampleSize, TrueRho, DF, MinSampleSizeForApprox,
        RhoExtreme, UseApprox,
        Sqrt2Pi
    end
    
    methods (Static)
     
        % References on Fisher's r to z:  Marascuilo1971, MengRosenthalRubin1992, SilverDunlap1987, ZimmermanZumboWilliams2003
        
        function z = Fisherrtoz(r)
            z = 0.5*log( (1+r) ./ (1 - r) );
        end
        
        function p = ApproxCDF(r,truerho,N)
            % Compute CDF(r) within rNoncentral having truerho & N.
            % truerho may be a vector
            p = zeros(size(truerho));
            zofr = rNoncentral.Fisherrtoz(r);
            se = 1 / sqrt(N-3);
            for iel=1:numel(truerho) 
                meanz = rNoncentral.Fisherrtoz(truerho(iel));
                p(iel) = normcdf( (zofr-meanz)/se );
            end
        end
        
        function r = Fisherztor(z)
            r = (exp(2*z)-1) / (exp(2*z)+1);
        end
        
        function r = ApproxInverseCDF(p,rho,N)
            z = norminv(p);
            meanz = rNoncentral.Fisherrtoz(rho);
            se = 1 / sqrt(N-3);
            r = rNoncentral.Fisherztor(z * se + meanz);
        end
        
    end
    
    methods
        
        function obj=rNoncentral(varargin)
            obj=obj@dContinuous('rNoncentral');
            obj.ParmTypes = 'ir';
            obj.DefaultParmCodes = 'fr';
            obj.NDistParms = 2;
            obj.CDFNearlyZero = 1e-6;
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            obj.RhoExtreme = 0.99;    % A bound to keep TrueRho from running away in searches.
            obj.MinSampleSizeForApprox = 120;
            obj.Sqrt2Pi = sqrt(2*pi);
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('rNoncentral:Constructor', ...
                        'rNoncentral constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.SampleSize = VerifyIntegerGE(obj,3,newparmvalues(1));
            obj.TrueRho = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newSampleSize = ifelse(ParmCodes(1)=='f', obj.SampleSize, obj.SampleSize+1);
            newnoncen = ifelse(ParmCodes(2)=='f', obj.TrueRho,1.1*obj.TrueRho);
            obj.ResetParms([newSampleSize newnoncen]);
        end
        
        function []=ReInit(obj)
            if ~((obj.TrueRho>=-obj.RhoExtreme)&&(obj.TrueRho<=obj.RhoExtreme))
                warning(['Resetting rNoncentral TrueRho ' num2str(obj.TrueRho) ' which must be between +/-' num2str(obj.RhoExtreme)]);
                if obj.TrueRho>obj.RhoExtreme
                    obj.TrueRho = obj.RhoExtreme;
                else
                    obj.TrueRho = -obj.RhoExtreme;
                end
            end
            obj.DF = obj.SampleSize - 2;
            obj.UseApprox = obj.SampleSize >= obj.MinSampleSizeForApprox;
            
            obj.Initialized = true;
            obj.LowerBound = -1 + eps;
            obj.UpperBound = 1 - eps;
            obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
            obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(3,Parms(1)) NumTrans.Bounded2Real(-obj.RhoExtreme,obj.RhoExtreme,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(3,Reals(1)) NumTrans.Real2Bounded(-obj.RhoExtreme,obj.RhoExtreme,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)  % Adapted from the function by Joshua Carmichael shown at the end of this file
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if obj.UseApprox
                thispdf = (obj.SampleSize-2) * (1-obj.TrueRho.^2)^((obj.SampleSize-1)/2) * (1-X.^2).^((obj.SampleSize-4)/2);
                thispdf = thispdf ./ (obj.Sqrt2Pi * (1-obj.TrueRho.*X).^(obj.SampleSize-3/2)) * obj.SampleSize.^(-1/2);
                thispdf = thispdf .* (1+ 1/4*(obj.TrueRho.*X+1)/(2*obj.SampleSize-1) + 9/16*(obj.TrueRho.*X+1).^2 / (2*obj.SampleSize-1)/(2*obj.SampleSize+1));
            else
                thispdf = (obj.SampleSize-2) * gamma(obj.SampleSize-1) * ((1-obj.TrueRho.^2).^((obj.SampleSize-1)/2)).* (1-X.^2).^((obj.SampleSize-4)/2);
                thispdf = thispdf ./ (obj.Sqrt2Pi * gamma(obj.SampleSize-1/2) * (1-obj.TrueRho.*X).^(obj.SampleSize-3/2));
                thispdf = thispdf .* (1+ 1/4*(obj.TrueRho.*X+1)/(2*obj.SampleSize-1) + 9/16*(obj.TrueRho.*X+1).^2 / (2*obj.SampleSize-1)/(2*obj.SampleSize+1));
            end
            thispdf(X>obj.UpperBound) = 0;
            thispdf(X<obj.LowerBound) = 0;
            %  thispdf(~isfinite(thispdf))     = 0;
        end
        
    end  % methods
    
end  % class rNoncentral

%   Here is this model for the PDF routine shown above.
%   function y = corrdist(r, ro, n)
%   % This function computes the probability density function for the
%   % correlation coefficient of a bivariate random variable.
%   %
%   % USAGES
%   % y = corrdist(r, ro, n)
%   %
%   % INPUT
%   % r:    Vector of possible correlation random variables, i.e. the values at
%   %       which the pdf is evaluated.
%   % ro:   The given (true) correlation coefficient, i.e. the population
%   %       correlation coefficient. length(ro) > 1 supported.
%   % n:    The number of samples in the correlated data. Only length(n) = 1
%   %       supported.
%   %
%   % OUTPUT
%   % y:    The probability density function for r, given ro, for n data
%   %       samples of a bivariate normal distribution.
%   %
%   %-----------------------------------------------------------------------
%   % Latest Edit: 11.June.2012
%   % Joshua D Carmichael
%   % josh.carmichael@gmail.com
%   %
%   % Original Author: Xu Cui, Stanford University (retrieved 11.June.2012)
%   %-----------------------------------------------------------------------
%
%   %accept vectorized inputs.
%   if(length(ro)> 1.5),
%       r   = repmat(r(:),1,length(ro));
%       ro  = repmat(ro(:)', length(r),1);
%   end;
%
%   if( n < 120 ),
%
%       y = (n-2) * gamma(n-1) * ((1-ro.^2).^((n-1)/2)).* (1-r.^2).^((n-4)/2);
%       y = y./ (sqrt(2*pi) * gamma(n-1/2) * (1-ro.*r).^(n-3/2));
%       y = y.* (1+ 1/4*(ro.*r+1)/(2*n-1) + 9/16*(ro.*r+1).^2 / (2*n-1)/(2*n+1));
%
%   else
%
%       y = (n-2) * (1-ro.^2)^((n-1)/2) * (1-r.^2).^((n-4)/2);
%       y = y./ (sqrt(2*pi) * (1-ro.*r).^(n-3/2)) * n.^(-1/2);
%       y = y.* (1+ 1/4*(ro.*r+1)/(2*n-1) + 9/16*(ro.*r+1).^2 / (2*n-1)/(2*n+1));
%
%   end;
%
%   y(r>1)              = 0;
%   y(r<-1)             = 0;
%   y(~isfinite(y))     = 0;
%
%   return;
