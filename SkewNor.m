classdef SkewNor < dContinuous
    % SkewNor(Location,Scale,Shape):  See the introduction at http:%azzalini.stat.unipd.it/SN/Intro/intro.html.
    
    properties(SetAccess = protected)
        Loc, Scale, Shape,
        Sqrt2OverPi, Delta, EX, VarX, Rho, Sqrt1minusRhoSqr,
        Standard_Normal, ZExtreme
        ComputingBounds  % Switch used to inhibit CDF integration when computing bounds
        InverseCDF02  % An arbitrary cutoff below which I integrate the PDF to get the CDF
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(0,Parms(2)) Parms(3)] ;
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(0,Reals(2)) Reals(3)];
        end
        
    end
    
    methods
        
        function obj=SkewNor(varargin)
            obj=obj@dContinuous('SkewNor');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.Standard_Normal = Normal(0,1);
            obj.Sqrt2OverPi = sqrt(2/pi);
            obj.ZExtreme = 30;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('SkewNor:Constructor', ...
                        'SkewNor constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Loc = newparmvalues(1);
            obj.Scale = newparmvalues(2);
            obj.Shape = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newloc   = ifelse(ParmCodes(1)=='f', obj.Loc, 1.02*obj.Loc);
            newscale = ifelse(ParmCodes(2)=='f', obj.Scale,0.98*obj.Scale);
            newshape = ifelse(obj.Shape==0,0.1,1.02*obj.Shape);
            newshape = ifelse(ParmCodes(3)=='f', obj.Shape,newshape);
            obj.ResetParms([newloc newscale newshape]);
        end
        
        function []=ReInit(obj)
            assert(obj.Scale>0,'SkewNor Scale must be positive.');
            obj.Delta = obj.Shape / sqrt(1 + obj.Shape^2);
            obj.EX = obj.Sqrt2OverPi*obj.Delta;
            obj.VarX = 1 - 2*obj.Delta^2/pi;
            obj.Rho = obj.Shape / sqrt( 1 + obj.Shape^2 );
            obj.Sqrt1minusRhoSqr = sqrt( 1 - obj.Rho^2 );
            obj.LowerBound = obj.Loc - obj.ZExtreme * obj.Scale;
            obj.UpperBound = obj.Loc + obj.ZExtreme * obj.Scale;
            obj.Initialized = true;
            obj.ComputingBounds = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            obj.InverseCDF02 = InverseCDF(obj,0.02);
            obj.ComputingBounds = false;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            Z = (X(InBounds) - obj.Loc) / obj.Scale;
            thispdf(InBounds) = 2 / obj.Scale * obj.Standard_Normal.PDF(Z) .* obj.Standard_Normal.CDF(obj.Shape*Z);
            %            for i=1:numel(X)
            %                if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %                    Z = (X(i) - obj.Loc) / obj.Scale;
            %                    thispdf(i) = 2 / obj.Scale * obj.Standard_Normal.PDF(Z) * obj.Standard_Normal.CDF(obj.Shape*Z);
            %                end
            %            end
        end
        
        function thiscdf=CDF(obj,X)
            % Using Owen's T function. See https://en.wikipedia.org/wiki/Skew_normal_distribution
            % This function produced negative CDFs eg
            % SkewNor(500.276,321.6111,25.0131).CDF(464) = -0.0019775
            % even though LowerBound = 345.77
            % In an attempt to fix this, I tried various cluges, including
            % the final cluge setting negative CDF values to zero.
            % But none of these cluges produced satisfactory results,
            % so I am falling back on the default integration CDF routine
            % for small X values, which seems more accurate though slower.
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            Xin = X(InBounds);
            Z = (X(InBounds) - obj.Loc) / obj.Scale;
            norcdf = normcdf(Z); %  obj.Standard_Normal.CDF(Z);
            wantpos = find(InBounds>0);
            for i=1:numel(Z)
                if ~obj.ComputingBounds && (Xin(i) < obj.InverseCDF02)
                    % Computing CDF by integration does not work until bounds have been set
                    thisone = CDF@dContinuous(obj,Xin(i));
                else
                    thisone = norcdf(i) - 2*TfnOwen(Z(i),obj.Shape);
                end
                thiscdf(wantpos(i)) = thisone;
            end
%             thiscdf(thiscdf<0) = 0;  % CLUGE HERE
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.Loc + obj.Scale*obj.Sqrt2OverPi*obj.Delta;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.Scale^2*(1 - 2*obj.Delta^2/pi);
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            CenMom3 = (4 - pi)/2 * obj.EX^3;  % / PowerOf(obj.VarX,1.5);
            if CenMom3 >= 0
                thisval = CenMom3^(1/3);
            else
                thisval = -(-CenMom3)^(1/3);
            end
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 3 + 2*(pi-3)*obj.EX^4/obj.VarX^2;
        end
        
        function thisval=Random(obj,varargin)
            % From http:%azzalini.stat.unipd.it/SN/faq-r.html, 2010-01-06:
            %Consider first the scalar SN case.
            %1. Sample u_0, u_1 having marginal distribution N(0,1) and correlation \obj.Rho.
            %   A simple way to achieve this is to generate u_0, v as independent N(0,1)
            %   variates and define u_1 = \obj.Rho \cdot u_0 + \sqrt(1-\obj.Rho^2)\cdot v.
            %2. Then
            %   z = (  u_1,  if u_0 > 0
            %       ( -u_1,  otherwise
            %   is a random number sampled from the SN distribution with obj.Shape parameter
            %   \alpha = \obj.Rho/\sqrt(1-\obj.Rho^2).
            %3. To change the location and obj.Scale from (0,1) to (a,b) with b>0, say, set y=a+b*z.
            %JOM notes: So we need to find \obj.Rho to produce the desired value of obj.Shape.
            %   obj.Shape = \obj.Rho/\sqrt(1-\obj.Rho^2)
            %   obj.Shape/\sqrt(1+obj.Shape^2) = \obj.Rho
            assert(obj.Initialized,UninitializedError(obj));
            thisval=zeros(varargin{:});
            for i=1:numel(thisval)
                u0 = obj.Standard_Normal.Random;
                v  = obj.Standard_Normal.Random;
                u1 = obj.Rho*u0 + obj.Sqrt1minusRhoSqr*v;
                if u0 > 0
                    RSN = u1;
                else
                    RSN = -u1;
                end
                thisval(i) = obj.Loc + obj.Scale*RSN;
            end
        end
        
        function parms = StartParmsMLE(obj,X)
            % Very seat-of-the-pants guesses based on some trial and error
            obs = prctile(X,[10 50 90]);
            thisMedian = obs(2);
            thissd = std(X);
            shapeguess = 10*(obs(1) + obs(3) - 2*thisMedian) / thissd;
            parms = [thisMedian 1.5*thissd shapeguess];
        end
        
    end  % methods
    
end  % class SkewNor

