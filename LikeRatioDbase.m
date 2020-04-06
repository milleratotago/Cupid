classdef LikeRatioDbase < dDiscrete
    % Parallel to LikeRatioCbase but for discrete data distributions,
    % but without any restrictions on H0Dist and H1Dist.
    % Note that there is no NBinsOrListOfX parameter, because LLR
    % can simply be computed for all of the discrete X's
    % By default, only parameters of the Data distribution are adjusted during estimation
    % (this makes sense to me).

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

        function obj=LikeRatioDbase(DataDist,H0Dist,H1Dist,sFamilyName)  % Constructor
            obj=obj@dDiscrete(sFamilyName);  % Inherited constructor
            if DataDist.DistType ~= 'd'
                error('Data distribution must be discrete to use LikeRatioD or LnLikeRatioD.');
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
            obj.DefaultParmCodes(obj.H0Start:obj.H1End) = 'f'; 
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
            % For each X, compute LLR
            Xs = obj.DataDist.DiscreteX;
            LLR = obj.LRorLnLR(Xs);
            [obj.DiscreteX, obj.DiscretePDF] = CollapseVals(LLR,obj.DataDist.DiscretePDF);
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.DiscreteCDF(end) = 1;
            obj.SetBinEdges;
            obj.StoredTablesInitialized = true;
            obj.NValues = numel(obj.DiscreteX);
            % Here I assume that DiscreteX is sorted.
            obj.LowerBound = obj.DiscreteX(1);
            obj.UpperBound = obj.DiscreteX(end);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
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

end  % LikeRatioDbase
