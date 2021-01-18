classdef ghHoag < dTransMono  % Lots of numerical problems with some parameter values--see unit tests.
    % ghHoag(A, B, g, h): From Hoaglin (1985) chapter, "Summarizing Shape Numerically: The g-and-h Distributions"
    % A & B>0 are location (median) and scale
    % g controls amount & direction of skew: g=0 is symmetric but leads to numerical errors so a cluge ensures abs(g)>eps
    %   typical values of g range between -1 and 1
    % h>0 controls heaviness of tails: h=0 is normal
    %   typical values of h range between 0 and 0.5
    % "A negative value of h is not impossible, but special treatment may be required because the monotonicity of Yh(z) fails for z^2 > -1/h."
    % "With g = 1 and h = 0, the distribution has the same shape as a lognormal distribution" Rousselet & Wilcox (2020)
    
    % Note: Unable to find reverse transformation for TransToPreTrans function:
    %   syms Z A B g h Trans
    %   transeqn = Trans == A + B * ( exp(g*Z) - 1 ) / g * exp(h*Z^2/2);
    %   solution = solve(transeqn,Z,'Real',true)
    %   > Warning: Unable to find explicit solution.
    % Thus, this distribution uses dTransMono.fzeroOpts.
    % I also tried to get an approximation for an improved starting point by dropping the minus 1:
    %     approx =  -(obj.g + (obj.g^2 + 2*obj.h*log(-(obj.g*(obj.A - TransX(i)))/obj.B))^(1/2))/obj.h;
    % alternative approximation:
    %     approx =  -(obj.g - (obj.g^2 + 2*obj.h*log(-(obj.g*(obj.A - TransX(i)))/obj.B))^(1/2))/obj.h;
    
    properties(Constant)
        minabsg = 0.0001; % min absolute value of g to avoid numerical errors if g too close to 0
        ming = -4;
        maxg =  4;
        minh =  0;
        maxh = 16;
    end
    
    properties(SetAccess = protected)
        A, B, g, h
    end
    
    methods
        
        function obj=ghHoag(A, B, g, h)
            zExtreme = 5;
            fixedNormal = Normal(0,1,zExtreme);
            fixedNormal.DefaultParmCodes = 'ff';
            obj=obj@dTransMono('ghHoag',fixedNormal);
            obj.fzeroOpts = optimset;
            obj.ResetParms([A B g h]);
            obj.AddParms(4,'rrrr');
            obj.PDFScaleFactorKnown = false;
            obj.ReInit;
        end

        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.A) ',' num2str(obj.B) ',' num2str(obj.g) ',' num2str(obj.h) ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.A = newparmvalues(end-3);
            obj.B = newparmvalues(end-2);
            obj.g = newparmvalues(end-1);
            if abs(obj.g) < obj.minabsg
                if obj.g < 0
                    obj.g = -obj.minabsg;
                else
                    obj.g = obj.minabsg;
                end
            end
            obj.h = newparmvalues(end);
            ClearBeforeResetParms(obj);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
%           obj.BasisRV.PerturbParms(ParmCodes);
            NewA = ifelse(ParmCodes(end-3)=='f', obj.A,obj.A*1.02);
            NewB = ifelse(ParmCodes(end-2)=='f', obj.B,obj.B*1.02);
            Newg = ifelse(ParmCodes(end-1)=='f', obj.g,obj.g*1.02);
            Newh = ifelse(ParmCodes(end  )=='f', obj.h,obj.h*1.02);
            obj.ResetParms([NewA NewB Newg Newh]);
        end
 
        function parmvals = TransParmValues(obj)
            parmvals = [obj.A obj.B obj.g obj.h];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            % Hoaglin (1985) p. 486, eqn 30b
            Y = ( exp(obj.g*PreTransX) - 1 ) / obj.g .* exp(obj.h*PreTransX.^2/2);
            TransX = obj.A + obj.B * Y;
            if TransX > realmax
                TransX = realmax;
            elseif TransX < -realmax
                TransX = -realmax;
            end
        end

        function PreTransX = TransToPreTrans(obj,TransX)
            if obj.UseSplineTransX
                PreTransX = spline(obj.SplineTransX,obj.SplineX,TransX);
            else
                PreTransX = zeros(size(TransX));
                for i=1:numel(TransX)
                    fn = @(x) (obj.PreTransToTrans(x) - TransX(i));
                    approx = 0;
%                   try
                        PreTransX(i) = fzero(fn,approx,obj.fzeroOpts);
%                   catch
%                      disp('fzero failed in TransToPreTrans');
%                      pause
%                   end % try
                end % for
            end % else
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = [Parms(end-3) Parms(end-2) NumTrans.Bounded2Real(obj.ming,obj.maxg,Parms(end-1)) NumTrans.Bounded2Real(obj.minh,obj.maxh,Parms(end))];
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = [Reals(end-3) Reals(end-2) NumTrans.Real2Bounded(obj.ming,obj.maxg,Reals(end-1)) NumTrans.Real2Bounded(obj.minh,obj.maxh,Reals(end))];
        end
        
        function x = Mean(obj)
            if (0 <= obj.h) && (obj.h <= 1)
                x = 1 / (obj.g * sqrt(1-obj.h)) * (exp(obj.g^2 / (2*(1-obj.h))) - 1);
                x = obj.A + obj.B * x;
            else
                x = Mean@dTransMono(obj);
            end
        end
        
        function x = RawMoment(obj,n)
            if (0 <= obj.h) && (obj.h <= 1/n) && (obj.A == 0)
                sum = 0;
                for i=0:n
                    term = (-1)^i * nchoosek(n,i) * exp( ((n-i)*obj.g)^2 / (2*(1-n*obj.h)) );
                    sum = sum + term;
                end
                x = 1 / (obj.g^n * sqrt(1-n*obj.h)) * sum;
                x = x * obj.B^n;
                % (exp(obj.g^2 / (2*(1-obj.h))) - 1);
            else
                x = RawMoment@dTransMono(obj,n);
            end
        end
        
    end  % methods
    
end  % class ghHoag



