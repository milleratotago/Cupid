classdef Rayleigh < dContinuous
    % Rayleigh(Sigma):  Sigma>0 is scale parameter
    
    properties(SetAccess = protected)
        Sigma, SigmaSqr
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = NumTrans.GT2Real(0,Parms(1));
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = NumTrans.Real2GT(0,Reals(1));
        end
        
        function thisval=RelSkewness()
            thisval = 0.631111;
        end
        
        function thisval=Kurtosis()
            thisval = 3.24509;
        end
        
    end
    
    methods
        
        function obj=Rayleigh(varargin)   % Constructor
            obj=obj@dContinuous('Rayleigh');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            
            % Handle constructor calls with different numbers of parameters.
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('Rayleigh:Constructor', ...
                        'Too many arguments passed to Rayleigh constructor.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Sigma = newparmvalues(1);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            if ~(ParmCodes(1)=='f')
                if obj.Sigma<5
                    Factor = 1.1;
                else
                    Factor = 0.9;
                end
                obj.ResetParms(Factor*obj.Sigma);
            end
        end
        
        function []=ReInit(obj)
            assert(obj.Sigma>=0,'Rayleigh Sigma must be >= 0.');
            obj.SigmaSqr = obj.Sigma^2;
            obj.Initialized = true;
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
            thispdf(InBounds) = X(InBounds) ./ ( exp(0.5*X(InBounds).^2/obj.SigmaSqr) * obj.SigmaSqr );
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = 1 - exp(-X(InBounds).^2*0.5/obj.SigmaSqr);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            thisval(InBounds) = obj.Sigma * sqrt( log( 1 ./ (1 - P(InBounds)).^2 ) );
            % for i=1:numel(P)
            %     thisval(i) = obj.Sigma * sqrt( log( 1 / (1 - P(i))^2 ) );
            % end
        end
        
        function thisval=Mean(obj)
            thisval = sqrt(pi/2) * obj.Sigma;
        end
        
        function thisval=Variance(obj)
            thisval = (2 - pi / 2) * obj.SigmaSqr;
        end
        
    end  % methods
    
end  % class Rayleigh





