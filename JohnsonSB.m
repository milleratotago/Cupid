classdef JohnsonSB < dContinuous
    % JohnsonSB(Location,Scale>0,Alpha1,Alpha2>0)
    % Warning: Parameter estimation does not work well (now skipped in utJohnsonSB).
    
    properties(SetAccess = protected)
        Location, Scale, Alpha1, Alpha2,
        MyZExtreme, Standard_Normal
    end
    
    methods
        
        function obj=JohnsonSB(varargin)
            obj=obj@dContinuous('JohnsonSB');
            obj.NDistParms = 4;
            obj.ParmTypes = 'rrrr';
            obj.DefaultParmCodes = 'rrrr';
            obj.Standard_Normal = Normal(0,1);
            obj.MyZExtreme = 10;  % May need adjustment depending on Alpha2
            switch nargin
                case 0
                case 4
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('JohnsonSB:Constructor', ...
                        'JohnsonSB constructor needs 0 or 4 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Location = newparmvalues(1);
            obj.Scale = newparmvalues(2);
            obj.Alpha1 = newparmvalues(3);
            obj.Alpha2 = newparmvalues(4);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newlocation = ifelse(ParmCodes(1)=='f', obj.Location,obj.Location+.1);
            newscale  = ifelse(ParmCodes(2)=='f', obj.Scale, 1.01*obj.Scale);
            newalpha1 = ifelse(ParmCodes(3)=='f', obj.Alpha1, 0.99*obj.Alpha1);
            newalpha2 = ifelse(ParmCodes(4)=='f', obj.Alpha2, 1.01*obj.Alpha2);
            obj.ResetParms([newlocation newscale newalpha1 newalpha2]);
        end
        
        function []=ReInit(obj)
            assert(obj.Scale>0,'JohnsonSB Scale must be > 0.');
            assert(obj.Alpha2>0,'JohnsonSB Alpha2 must be > 0.');
            obj.LowerBound = obj.Location + obj.Scale / ( 1 + exp( -(-obj.MyZExtreme-obj.Alpha1)/obj.Alpha2 ) );
            obj.UpperBound = obj.Location + obj.Scale / ( 1 + exp( -( obj.MyZExtreme-obj.Alpha1)/obj.Alpha2 ) );
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2)) Parms(3) NumTrans.GT2Real(eps,Parms(4))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2)) Reals(3) NumTrans.Real2GT(eps,Reals(4))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            X(InBounds) = (X(InBounds) - obj.Location) / obj.Scale;
            Z = obj.Alpha1 + obj.Alpha2 * log(X(InBounds) ./ (1 - X(InBounds)));
            NorPDF = obj.Standard_Normal.PDF(Z);
            thispdf(InBounds) = obj.Alpha2 ./ (obj.Scale*X(InBounds).*(1-X(InBounds))) .* NorPDF;
            %            for i=1:numel(X)
            %                if (X(i)>=obj.LowerBound) && (X(i)<=obj.UpperBound)
            %                    X(i) = (X(i) - obj.Location) / obj.Scale;
            %                    Z = obj.Alpha1 + obj.Alpha2 * log(X(i) / (1 - X(i)));
            %                    NorPDF = obj.Standard_Normal.PDF(Z);
            %                    thispdf(i) = obj.Alpha2 / (obj.Scale*X(i)*(1-X(i))) * NorPDF;
            %                end
            %            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            X(InBounds) = (X(InBounds) - obj.Location) / obj.Scale;
            XOverOneMinusX = X(InBounds) ./ (1 - X(InBounds));
            Z = obj.Alpha1 + obj.Alpha2 * log(XOverOneMinusX);
            thiscdf(InBounds) = CDF(obj.Standard_Normal,Z);
            %             for i=1:numel(X)
            %                 if X(i)<=obj.LowerBound
            %                 elseif X(i)>=obj.UpperBound
            %                     thiscdf(i) = 1;
            %                 else
            %                     X(i) = (X(i) - obj.Location) / obj.Scale;
            %                     XOverOneMinusX = X(i) / (1 - X(i));
            %                     Z = obj.Alpha1 + obj.Alpha2 * log(XOverOneMinusX);
            %                     thiscdf(i) = CDF(obj.Standard_Normal,Z);
            %                 end
            %             end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            Z = Random(obj.Standard_Normal,varargin{:});
            thisval = obj.Location + obj.Scale ./ ( 1 + exp( -(Z-obj.Alpha1)/obj.Alpha2 ) );
        end
        
    end  % methods
    
end  % class JohnsonSB



