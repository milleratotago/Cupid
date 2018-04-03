classdef RNGamma < dContinuous
    % RNGamma(N,Rate) where N is a real number.
    
    properties(SetAccess = private)    % These properties can only be set by the methods of this class.
    end
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        N, Rate, LnGammaOfN, Standard_Exponential
    end
    
    methods
        
        function obj=RNGamma(varargin)   % Constructor
            obj=obj@dContinuous('RNGamma');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.Standard_Exponential = Exponential(1);
            obj.IntegralPDFXmuNAbsTol = 10*obj.IntegralPDFXmuNAbsTol;  % For integrating PDF
            obj.IntegralPDFXmuNRelTol = 10*obj.IntegralPDFXmuNRelTol;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('RNGamma:Constructor', ...
                        'RNGamma constructor must receive 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.N = newparmvalues(1);
            obj.Rate = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, increase parmN and decrease parmRate
            newN    = ifelse(ParmCodes(1)=='f', obj.N,    1.1*obj.N);
            newRate = ifelse(ParmCodes(2)=='f', obj.Rate, 0.9*obj.Rate);
            obj.ResetParms([newN newRate]);
        end
        
        function []=ReInit(obj)
            assert(obj.N>0,'RNGamma N must be > 0.');
            assert(obj.Rate>0,'RNGamma Rate must be > 0.');
            obj.LnGammaOfN = gammaln(obj.N);
            obj.LowerBound = obj.XNearlyZero;
            obj.UpperBound = obj.Standard_Exponential.InverseCDF(obj.CDFNearlyOne) / obj.Rate * obj.N;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(0,Parms(1)) NumTrans.GT2Real(0,Parms(2))] ;
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(0,Reals(1)) NumTrans.Real2GT(0,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = gampdf(X(InBounds),obj.N,1/obj.Rate);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = gamcdf(X(InBounds),obj.N,1/obj.Rate);
        end
        
        
        % I think the next two only work for integer obj.N:
        %
        %         function thisval=MGF(obj,Theta)
        %             assert(obj.Initialized,UninitializedError(obj));
        %             Power = (1 - Theta * obj.Rate)^obj.N;
        %             thisval = 1 / Power;
        %         end
        %
        %         function thisval=RawMoment(obj,I)
        %             assert(obj.Initialized,UninitializedError(obj));
        %             Sum1 = 0;
        %             Sum2 = 0;
        %             for J = 2:obj.N + I - 1
        %                 Add = log(J);
        %                 Sum1 = Sum1 + Add;
        %                 if J <= obj.N-1
        %                     Sum2 = Sum2 + Add;
        %                 end
        %             end
        %             thisval = exp(Sum1 - Sum2 - I * log(obj.Rate) );
        %         end
        %
        %
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            %
            % Tested successfully with N=2.5, Rate = 3.33;
            % Not tested with N < 1.
            %
            % Information from Devroye 1986:
            % Density Gamma(N,Rate).
            %  N>0 is shape parameter.
            %  Rate>0 is scale parameter.
            % if X is Gamma(N,Rate) then cX is Gamma(N,bc) for c>0.
            %
            % Gamma(N,1) is abbreviated Gamma(N).
            %
            % This generator produces Gamma(N), and then Divides by Rate at end.
            %
            % There are two methods here, both from Devroye.
            %
            % With N>=1, use the algorithm of Best, 1978, page 410.
            % With 0<N<1, use modified (apparently by Devroye) version of algorithm
            %   of Vaduva (1977), page 415.  Be careful of typo in book:
            %   on page 415,     Until Z + E <= d + X  should be >= instead of <=
            %     bb, bc, bU, bV, bW, bX, bY, bZ,    % Used by Best
            %     vc, vd, vZ, vE, vX,                % Used by (modified) Vaduva
            thisval = zeros(varargin{:});
            for iel=1:numel(thisval)
                if obj.N == 1
                    GammaA1 = obj.Standard_Exponential.Random;
                elseif obj.N > 1    % Best's algorithm to generate Gamma(N,1)
                    bb = obj.N - 1;
                    bc = 3*obj.N - 0.75;
                    Accept = false;
                    while ~Accept
                        bU = rand;
                        bV = rand;
                        bW = bU * (1 - bU);
                        if bW == 0
                            Accept = false;
                        else
                            bY = sqrt(bc/bW) * (bU - 0.5);
                            bX = bb + bY;
                            if bX >= 0
                                bZ = 64 * bW * bW^2 * bV^2;
                                % The following is the acceptance criterion in the original Devroye,
                                % but it fails when bb=0. The corrected version (see Erratum) is below.
                                % The correction is actually irrelevant to this code because I use the
                                % exponential generator when N=1:    "if N = 1 Then GammaA1 = Standard_Exponential.Random"
                                %          Accept = ( bZ <= 1 - 2 * Sqr(bY) / bX ) or
                                %                    ( log(bZ) <= 2 * (bb * log(bX / bb) - bY) );
                                Accept = ( bZ <= 1 - 2 * bY^2 / bX ) || ...
                                    ( (bb > 0) && ( log(bZ) <= 2 * (bb * log(bX / bb) - bY) ) ) || ...
                                    ( (bb == 0) && ( log(bZ) <= -2 * bY) );
                            end
                        end
                    end  % while ~Accept
                    GammaA1 = bX;
                    % Best's algorithm
                else    % Devroye's modification of Vaduva (1977) algorithm
                    vc = 1.0 / obj.N;
                    vd = obj.N^(obj.N/(1-obj.N)) * (1-obj.N);
                    Accept = false;
                    while ~Accept
                        vZ = obj.Standard_Exponential.Random;
                        vE = obj.Standard_Exponential.Random;
                        vX = vZ^vc;  % vX is Weibull(N)
                        Accept = vZ + vE >= vd + vX;
                    end
                    GammaA1 = vX;
                    %  Writeln('Use Vaduva');
                end % modified Vaduva
                thisval(iel) = GammaA1 / obj.Rate; % Return scaled Gamma(N,1)
            end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.N / obj.Rate;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.N / obj.Rate^2;
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 2 / sqrt(obj.N);
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 3 + 6.0 / obj.N;
        end
        
    end  % methods
    
end  % class RNGamma


