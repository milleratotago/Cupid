classdef ghHoag < dTransMono  % Lots of numerical problems with some parameter values--see unit tests.
    % ghHoag(A, B, g, h): From Hoaglin (1985) chapter, "Summarizing Shape Numerically: The g-and-h Distributions"
    % A & B>0 are location (median) and scale
    % g>0 controls amount & direction of skew: g=0 is symmetric but leads to numerical errors so a cluge ensures abs(g)>eps
    % h>0 controls heaviness of tails: h=0 is normal
    % "A negative value of h is not impossible, but special treatment may be required because the monotonicity of Yh(z) fails for z^2 > -1/h."
    
    properties(SetAccess = protected)
        A, B, g, h
        ming  % min absolute value of g to avoid numerical errors if g==0
        maxg
        minh, maxh
    end
    
    methods
        
        function obj=ghHoag(A, B, g, h)
            zExtreme = 5;
            fixedNormal = Normal(0,1,zExtreme);
            fixedNormal.DefaultParmCodes = 'ff';
            obj=obj@dTransMono('ghHoag',fixedNormal);
            obj.ming = 0.00001;
            obj.maxg = 2;
            obj.minh = -1;
            obj.maxh = 1;
            obj.ResetParms([A B g h]);
            obj.AddParms(4,'rrrr');
            obj.PDFScaleFactorKnown = false;
            obj.ReInit;
        end

        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.A) ',' num2str(obj.B) ',' num2str(obj.g) ',' num2str(obj.h) ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.A = newparmvalues(1);
            obj.B = newparmvalues(2);
            obj.g = newparmvalues(3);
            if abs(obj.g) < obj.ming
                if obj.g < 0
                    obj.g = -obj.ming;
                else
                    obj.g = obj.ming;
                end
            end
            obj.h = newparmvalues(4);
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
            % p. 486, eqn 30b
            Y = (exp(obj.g*PreTransX)-1) / obj.g .* exp(obj.h*PreTransX.^2/2);
            TransX = obj.A + obj.B * Y;
            if TransX > realmax
                TransX = realmax;
            elseif TransX < -realmax
                TransX = -realmax;
            end
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



