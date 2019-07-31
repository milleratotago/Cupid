classdef LikeRatioCbase < dContinuous
    % A base class used by both LikeRatioC and LnLikeRatioC
    % The basic idea is to use a spline approximation for the distribution of (ln) likelihood ratio values
    % that would be obtained when testing between two hypothesized distributions, H0Dist and H1Dist.
    % The data come from a distribution DataDist.
    % Restrictions:
    %    Larger values of the basis distribution must always relatively favor H1 over H0.
    %    Likelihood ratio values should be monotonic in the X values of the data distribution.
    %    There are problems when the data distribution includes values that are impossible under H0 or H1.
    % Note: Spline fails for cases where there are not enough unique values, e.g. where H0 & H1 are uniform so
    % that the ratio of pdfs is constant.

    properties(SetAccess = protected)
        DataDist
        H0Dist
        H1Dist

        % Indices of parameters & parmcodes for each of the three basis distributions
        % within the unidimensional vectors of parameters and parmcodes for the overall distribution.
        DataStart, DataEnd
        H0Start, H0End
        H1Start, H1End
    end

    methods(Abstract)
        thisval = LRorLnLR(obj,X)  % Descendants define this function to obtain the distribution
                                   % of likelihood scores, -2*Ln(likelihood) scores.
                                   % LRorLnLR must be vectorized.
    end  % Abstract methods

    methods

        function obj=LikeRatioCbase(DataDist,H0Dist,H1Dist,NBinsOrListOfX,sFamilyName)  % Constructor
            obj=obj@dContinuous(sFamilyName);  % Inherited constructor
            if DataDist.DistType ~= 'c'
                error('Data distribution must be continuous to use LikeRatioC or LnLikeRatioC.');
            end
            obj.DataDist = DataDist;
            obj.H0Dist = H0Dist;
            obj.H1Dist = H1Dist;
            obj.DataStart = 1;
            obj.DataEnd = obj.DataDist.NDistParms;
            obj.H0Start = obj.DataEnd + 1;
            obj.H0End = obj.H0Start + obj.H0Dist.NDistParms - 1;
            obj.H1Start = obj.H0End + 1;
            obj.H1End = obj.H1Start + obj.H1Dist.NDistParms - 1;
            obj.NDistParms = obj.DataDist.NDistParms + obj.H0Dist.NDistParms + obj.H1Dist.NDistParms;
            obj.DefaultParmCodes = [obj.DataDist.DefaultParmCodes obj.H0Dist.DefaultParmCodes obj.H1Dist.DefaultParmCodes];
            obj.SplineCDFsXSpec = NBinsOrListOfX;
            obj.ReInit;
        end

        function PerturbParms(obj,ParmCodes)
            obj.DataDist.PerturbParms(ParmCodes(obj.DataStart:obj.DataEnd));
            obj.H0Dist.PerturbParms(ParmCodes(obj.H0Start:obj.H0End));
            obj.H1Dist.PerturbParms(ParmCodes(obj.H1Start:obj.H1End));
            obj.ResetParms([obj.DataDist.ParmValues obj.H0Dist.ParmValues obj.H1Dist.ParmValues]);
        end
        
        function BuildMyName(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            obj.StringName = [obj.FamilyName '(' obj.DataDist.StringName ',' obj.H0Dist.StringName ',' obj.H1Dist.StringName ')'];
        end
        
        function ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj);
            obj.DataDist.ResetParms(newparmvalues(obj.DataStart:obj.DataEnd));
            obj.H0Dist.ResetParms(newparmvalues(obj.H0Start:obj.H0End));
            obj.H1Dist.ResetParms(newparmvalues(obj.H1Start:obj.H1End));
            obj.ReInit;
        end
        
        function [] = ReInit(obj)
            % Vary X over the range of DataDist
            % For each X, compute LLR and CDFx
            
            % Compute the X's at which to determine LikeRatio values:
            if numel(obj.SplineCDFsXSpec) == 1
                StepSize = (obj.DataDist.UpperBound - obj.DataDist.LowerBound) / obj.SplinePDFsXSpec;
                Xs = obj.DataDist.LowerBound+StepSize:StepSize:obj.DataDist.UpperBound-StepSize;
                % Look for bounds on the Xs such that LLR is a real number.
                WorkingLowerBound = BoundSearch(0,0.0001);
                WorkingUpperBound = BoundSearch(1,-0.0001);
                Xs = [WorkingLowerBound Xs WorkingUpperBound];
            else
                Xs = obj.SplineCDFsXSpec;
            end
            
            % Compute the (Ln)LikeRatio values at those Xs and select out unique ones:
            LLR = obj.LRorLnLR(Xs);
            FailAt = isinf(LLR) | isnan(LLR);
            Xs = Xs(~FailAt);
            LLR = LLR(~FailAt);
            [uniqueLLR, indices] = unique(LLR);
            correspondingXs = Xs(indices);
            obj.SplineCDFsXs = uniqueLLR;
            obj.SplineCDFs = obj.DataDist.CDF(correspondingXs);
           
            obj.LowerBound = min(obj.SplineCDFsXs);
            obj.UpperBound = max(obj.SplineCDFsXs);
            obj.CDFSplineInfo = spline(uniqueLLR,obj.SplineCDFs);
            obj.UseSplineCDF = true;
            obj.HaveSplineCDFs = true;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end

            function WorkingBound = BoundSearch(Startp,Increment)
                Found = false;
                while ~Found
                   Startp = Startp + Increment;
                   WorkingBound = obj.DataDist.InverseCDF(Startp);
                   thisLLR = obj.LRorLnLR(WorkingBound);
                   Found = ~ (isnan(thisLLR) || isinf(thisLLR));
                end
            end

        end

        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.DataDist.ParmsToReals(Parms(obj.DataStart:obj.DataEnd)) obj.H0Dist.ParmsToReals(Parms(obj.H0Start:obj.H0End)) obj.H1Dist.ParmsToReals(Parms(obj.H1Start:obj.H1End))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.DataDist.RealsToParms(Reals(obj.DataStart:obj.DataEnd)) obj.H0Dist.RealsToParms(Reals(obj.H0Start:obj.H0End)) obj.H1Dist.RealsToParms(Reals(obj.H1Start:obj.H1End))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.DataDist.ParmValues obj.H0Dist.ParmValues obj.H1Dist.ParmValues];
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            randXs = obj.DataDist.Random(varargin{:});
            thisval = obj.LRorLnLR(randXs);
        end
        
    end  % methods

end  % classdef
