classdef Geary < dContinuous
    % Geary(SampleSize)
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        SampleSize
        SqrtSampleSize
        Standard_Normal
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = NumTrans.GT2Real(1,Parms(1));
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = NumTrans.Real2GT(1,Reals(1));
        end
        
        function thisA = ObservedA(X)
            % Compute the observed value of Geary's A for the data in X.
            X = X(:);  % Make sure X is a vector.
            N = numel(X);
            XBar = mean(X);
            Dev = X - XBar;
            SumAbsDev = sum(abs(Dev));
            SumSqrDev = sum(Dev.^2);
            thisA = SumAbsDev / N / sqrt(SumSqrDev/N);
        end
        
        function ManyAs = SimulatedAs(SampleSize, NSims)
            % Generate many values of Geary's A by simulation.
            ManyAs = zeros(1,NSims);
            for i=1:NSims
                x = randn(SampleSize,1);
                ManyAs(i) = Geary.ObservedA(x);
            end
        end
        
    end
    
    methods
        
        function obj=Geary(SampleSize)
            obj=obj@dContinuous('Geary');
            obj.ParmTypes = 'I';
            obj.DefaultParmCodes = 'I';
            obj.NDistParms = 1;
            obj.Standard_Normal = Normal(0,1);
            obj.Initialized = true;
            obj.ResetParms(SampleSize);
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.SampleSize = newparmvalues(1);
            obj.SqrtSampleSize = sqrt(obj.SampleSize);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            New = ifelse(ParmCodes(1)=='f', obj.SampleSize, obj.SampleSize+1);
            obj.ResetParms(New);
        end
        
        function []=ReInit(obj)
            obj.LowerBound = 0;
            obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thiscdf = CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            ZofGeary = (X(InBounds) - 0.7979) * obj.SqrtSampleSize / 0.2123;
            thiscdf(InBounds) = obj.Standard_Normal.CDF(ZofGeary);
        end
        
        function thisval = InverseCDF(obj,P)
            ZofP = obj.Standard_Normal.InverseCDF(P);
            thisval = ZofP * 0.2123 / obj.SqrtSampleSize + 0.7979;
        end
        
    end  % methods
    
end  % class Geary

