classdef r < dContinuous
    
    properties(SetAccess = protected)
        SampleSize, DF, GamHalf,  Gam1, Gam2, Parentt
    end
    
    methods
        
        function obj=r(varargin)   % Constructor
            obj=obj@dContinuous('r');
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            obj.NDistParms = 1;
            obj.GamHalf = gamma(1/2);
            obj.Parentt = t;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('r:Constructor', ...
                        'r constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.SampleSize = VerifyIntegerGE(obj,3,newparmvalues(1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            if ~(ParmCodes(1)=='f')
                obj.ResetParms(obj.SampleSize+1);
            end
        end
        
        function []=ReInit(obj)
            obj.DF = obj.SampleSize - 2;
            obj.Gam1 = gamma( (obj.SampleSize - 1) / 2 );
            obj.Gam2 = gamma( (obj.SampleSize - 2) / 2 );
            ResetParms(obj.Parentt,obj.DF);
            obj.LowerBound = -1;
            obj.UpperBound = 1;
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = -obj.LowerBound;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.GT2Real(3,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(3,Reals(1));
        end
        
        function thispdf=PDF(obj,X)    % JKB
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) =  obj.Gam1 ./ (obj.GamHalf*obj.Gam2) * (1-X(InBounds).^2).^((obj.SampleSize-4)/2);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            Equivt = X(InBounds) .* sqrt( obj.DF ./ (1 - X(InBounds).^2) );
            thiscdf(InBounds) = obj.Parentt.CDF(Equivt);
        end
        
        function thisval=Mean(obj)
            thisval = 0;
        end
        
        function thisval=rRV.InverseCDF(P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            % Ulrich & Miller p-curve paper says that the inverse cdf is:
            % \FoInv(p) = \tanh \left[ \Phi^{-1}(p) \cdot \sqrt{n-3}  \right].
            thisval = tanh( norminv(p) * sqrt(n-3)  );
%             for iel=1:numel(P)
%                 if InBounds(iel)
%                     Equivt = obj.Parentt.InverseCDF(P(iel));
%                     Neg = Equivt < 0;
%                     r = sqrt( Equivt^2/obj.DF/(1+Equivt^2/obj.DF) );
%                     if Neg
%                         r = -r;
%                     end
%                     thisval(iel) = r;
%                 end
%             end
        end
        
    end  % methods
    
end  % class r

