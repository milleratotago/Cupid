classdef clOutlierModel < handle
    % Model scores in N separate conditions "contaminated" by a common outlier process
    % that "contaminates" scores in one of NTypes possible ways (see Enum below).
    % In all cases, the parameters of an object of this class are
    % stored in the order: TrueDists{i} parms, prOutlier, Contam parms.
    
    properties(Constant)
        NTypes  = 3;  % N of types enumerated here.
        % Enum of different types of contamination:
        Shift   = 1;  % Contamination RV adds to true score
        Stretch = 2;  % Contamination RV multiplies times true score
        Replace = 3;  % Contamination RV replaces true score
        
        TypeLabel = {'Shift', 'Stretch', 'Replace'};  % Labels used for file naming, etc.
        
    end
    
    properties
        
        % "Input" properties defining the model:
        TrueDists  % Cell array of Cupid distributions of true ("uncontaminated") scores in the N conditions.
        prOutlier  % Outlier probability, constant across conditions.
        Contam     % Cupid distribution representing trial-to-trial variation in the extent of contamination.
        ConType    % Type of contamination: Shift, Stretch, or Replace
        
        % Optional input parameters:
        
        % Properties derived from the input properties:
        NConds     % N of conditions described by the model
        ObsDists   % Cell array of Cupid distributions of observed scores in the N conditions, including outliers.
        
        DefaultParmCodes
        SearchOptions
        SplinePDFs
        SplineCDFs
        
        % Indices of parameters within the overall parameter list:
        TrueParmIndices   % Cell array 1:NConds with indices of each distribution's parameters
        prOutlierIndex    % Index of the prOutlier parameter
        ContamParmIndices % Vector of indices of the parameters for the Contam distribution
        TotalNParms
        
    end
    
    properties (SetAccess = public)
       minPrOutlier, maxPrOutlier % Boundaries on prOutlier used in parameter searches.
    end

    methods
        
        function obj = clOutlierModel(TrueDists,prOutlier,Contam,ConType,varargin)
            % Create object with specified values of input parameters.
            % TrueDists is either a single Cupid distribution or a cell array of
            %   Cupid distributions (1 per condition).
            assert(ConType >= 1 && ConType <= clOutlierModel.NTypes && ConType == floor(ConType),'Illegal ConType');
            assert(prOutlier > 0 && prOutlier < 1,'prOutlier must be between 0 & 1');
            [obj.SearchOptions, varargin] = ExtractNameVali('SearchOptions',optimset,varargin);
            [obj.SplinePDFs, varargin] = ExtractNameVali({'SplinePDF','SplinePDFs'},0,varargin);
            [obj.SplineCDFs, varargin] = ExtractNameVali({'SplineCDF','SplineCDFs'},0,varargin);
            assert(numel(varargin)==0,'Unrecognized argument');
            if iscell(TrueDists)
                obj.TrueDists = TrueDists;
            else
                % Allow user to specify a single true distribution simply as a Cupid
                % distribution rather than as a 1-position cell array.
                obj.TrueDists = {TrueDists};
            end
            obj.prOutlier = prOutlier;
            obj.minPrOutlier = min(0.00001,prOutlier);  % Defaults can be overridden
            obj.maxPrOutlier = max(0.20,prOutlier);
            obj.Contam = Contam;
            obj.ConType = ConType;
            obj.NConds = numel(TrueDists);
            obj.MakeObsDists;
            
            % Construct default parmcodes
            sParmCodes = '';
            for iCond = 1:obj.NConds
                sParmCodes = [sParmCodes obj.TrueDists{iCond}.DefaultParmCodes]; %#ok<AGROW>
            end
            sParmCodes = [sParmCodes 'r' obj.Contam.DefaultParmCodes];
            obj.DefaultParmCodes = sParmCodes';  % ' to match with column vector whose size is checked in fminsearcharb
            
            % Compute the indices of all parameters within the overall parameter lists
            UsedParms = 0;
            obj.TrueParmIndices = cell(obj.NConds,1);
            for iCond=1:obj.NConds
                obj.TrueParmIndices{iCond} = (UsedParms+1):(UsedParms+obj.TrueDists{iCond}.NDistParms);
                UsedParms = UsedParms + obj.TrueDists{iCond}.NDistParms;
            end
            obj.prOutlierIndex  = UsedParms + 1;
            obj.ContamParmIndices = (UsedParms+2):UsedParms+1+obj.Contam.NDistParms;
            obj.TotalNParms = obj.ContamParmIndices(end);
            
        end
        
        function [] = MakeOneObsDist(obj,iCond)
            % Initialize the observed distribution for condition iCond the current parameter values.
            switch obj.ConType
                case obj.Shift
                    obj.ObsDists{iCond} = ContamShift(obj.TrueDists{iCond},obj.prOutlier,obj.Contam);
                case obj.Stretch
                    obj.ObsDists{iCond} = ContamStretch(obj.TrueDists{iCond},obj.prOutlier,obj.Contam);
                case obj.Replace
                    obj.ObsDists{iCond} = Mixture(1-obj.prOutlier,obj.TrueDists{iCond},obj.Contam);
            end
            if obj.SplinePDFs > 0
                obj.ObsDists{iCond}.UseSplinePDFOn(obj.SplinePDFs);
            end
            if obj.SplineCDFs > 0
                obj.ObsDists{iCond}.UseSplineCDFOn(obj.SplineCDFs);
            end
        end
        
        function [] = MakeObsDists(obj)
            % Initialize the observed distributions with the current parameter values.
            obj.ObsDists = cell(obj.NConds,1);
            for iCond=1:obj.NConds
                obj.MakeOneObsDist(iCond);
            end
        end
        
        function Observed = Random(obj,varargin)
            % Generate random observed scores with current model parameters.
            % Observed is a cell array, 1:NConds, of scores in the different conditions.
            Observed = cell(obj.NConds,1);
            for iCond = 1:obj.NConds
                Observed{iCond} = obj.ObsDists{iCond}.Random(varargin{:});
            end
        end
        
        function ll = LnLikelihood(obj,X)
            % Compute likelihood of the data in X with current model parameters.
            % X is a cell array, 1:NConds, of scores in the different conditions.
            perCondll = zeros(obj.NConds,1);
            for iCond = 1:obj.NConds
                perCondll(iCond) = obj.ObsDists{iCond}.LnLikelihood(X{iCond});
            end
            ll = sum(perCondll);
        end
        
        function [EndingVals,fval,exitflag,output]=EstML(obj,X,varargin)
            % Estimate distribution parameters by maximum likelihood [i.e., minimize -log(likelihood)]
            % X is a cell array, 1:NConds, of observed scores in the different conditions.
            % exitflags: 1 converged, 0 N of iterations exceeded, -1 terminated by output function.
            if numel(varargin)<1
                ParmCodes = obj.DefaultParmCodes;
            else
                ParmCodes = varargin{1};
            end
            holdWarnStates = true(obj.NConds,1);
            for iCond=1:obj.NConds
               holdWarnStates(iCond) = obj.ObsDists{iCond}.SkipImpossibleWarn;
               obj.ObsDists{iCond}.SkipImpossibleWarn = true;
            end
            Real2ParmFn = @obj.RealsToParms;
            Parm2ReamFn = @obj.ParmsToReals;
            ErrFn = @MyErrFunc;
            StartingVals = ParmValues(obj);
            [EndingVals,fval,exitflag,output] = fminsearcharb(ErrFn,StartingVals,Real2ParmFn,Parm2ReamFn,ParmCodes,obj.SearchOptions);
            for iCond=1:obj.NConds
               obj.ObsDists{iCond}.SkipImpossibleWarn = holdWarnStates(iCond);
            end
            obj.ResetParms(EndingVals);
            
            function thiserrval=MyErrFunc(realParms)
                ResetParms(obj,realParms)
                thiserrval = -LnLikelihood(obj,X);
            end
            
        end
        
        function [EndingVals,fval,exitflag,output]=EstMLiter(obj,X)
            % Estimate distribution parameters by maximum likelihood [i.e., minimize -log(likelihood)]
            % X is a cell array, 1:NConds, of observed scores in the different conditions.
            % exitflags: 1 converged, 0 N of iterations exceeded, -1 terminated by output function.
            % *** ParmCodes IS NOT YET SUPPORTED
            % This version uses an iterative procedure that starts with a given PrOutlier and Contam parms:
            %   adjust the TrueParms for each distribution separately using the current PrOutlier and Contam parms
            %   adjust PrOutlier and Contam parms using the current TrueParms for each distribution
            %   check for convergence
            
            fvalConvergenceCutoff = 0.01;
            parmConvergenceCutoff = 0.01;
            MaxIter = 200;
            
            NContamParms = obj.Contam.NDistParms;
            FixedContamParmCodes = repmat('f',NContamParms,1);
            FreeContamParmCodes = repmat('r',NContamParms+1,1);   % extra 1 is prOutlier
            NTrueParms = numel(obj.DefaultParmCodes)-NContamParms-1;  % extra 1 is prOutlier
            FixedDistParmCodes = repmat('f',NTrueParms,1);
            Converged = false;
            iter = 0;
            while (iter<MaxIter) & ~Converged
                iter = iter + 1
                % Adjust the TrueParms for each distribution separately using the current PrOutlier and Contam parms.
                % Remember that the ObsDists are mixtures of true and contam dists & there is a mixture prob parameter
                ests = [];
                for iCond=1:obj.NConds
                    obj.MakeOneObsDist(iCond);
                    DistParmCodes = obj.ObsDists{iCond}.DefaultParmCodes;
                    DistParmCodes(end-NContamParms+1:end) = FixedContamParmCodes;
                    [~,est1] = obj.ObsDists{iCond}.EstML(X{iCond},DistParmCodes);
                    trueEsts = est1(DistParmCodes=='r');
                    ests = [ests trueEsts]; %#ok<AGROW>
                end % for iCond
%                 obj.ResetParms([ests obj.prOutlier obj.Contam.ParmValues]);
                
                % Adjust PrOutlier and Contam parms using the current TrueParms for each distribution.
                % First parameter (free) is mixture probability.
                [EndingVals,fval] = obj.EstML(X,[FixedDistParmCodes; FreeContamParmCodes]);
                [EndingVals' fval]
%                 obj.ResetParms(EndingVals);
                if iter > 1
                    parmChange = norm(EndingVals-prevEndingVals);
                    fvalChange = prevfval - fval;
                    Converged = (fvalChange < fvalConvergenceCutoff) && (parmChange < parmConvergenceCutoff);
                end
                prevEndingVals = EndingVals;
                prevfval = fval;
                
            end % while ~Converged
            exitflag = Converged;
            output = 'unknown';
        end
        
        
        function ResetParms(obj,newparmvalues)
            for iCond=1:obj.NConds
                theseIndices = obj.TrueParmIndices{iCond};
                obj.TrueDists{iCond}.ResetParms(newparmvalues(theseIndices));
            end
            obj.prOutlier = newparmvalues(obj.prOutlierIndex);
            obj.Contam.ResetParms(newparmvalues(obj.ContamParmIndices));
            obj.MakeObsDists;
        end
        
        function parms = ParmValues(obj)
            % Retrieve the current parameter values
            parms = zeros(obj.TotalNParms,1);
            for iCond=1:obj.NConds
                theseIndices = obj.TrueParmIndices{iCond};
                parms(theseIndices) = obj.TrueDists{iCond}.ParmValues;
            end
            parms(obj.prOutlierIndex) = obj.prOutlier;
            parms(obj.ContamParmIndices) = obj.Contam.ParmValues;
            
        end
        
        function realList = ParmsToReals(obj,Parms,~)
            % Return a list of the model parameter values transformed into reals as appropriate for fminsearch.
            % These come out in the order: TrueDists parms, prOutlier, Contam parms.
            realList = zeros(obj.TotalNParms,1);
            for iCond = 1:obj.NConds
                theseIndices = obj.TrueParmIndices{iCond};
                realList(theseIndices) = obj.TrueDists{iCond}.ParmsToReals(Parms(theseIndices));
            end
            realList(obj.prOutlierIndex) = NumTrans.Bounded2Real(obj.minPrOutlier,obj.maxPrOutlier,Parms(obj.prOutlierIndex));
            realList(obj.ContamParmIndices) = obj.Contam.ParmsToReals(Parms(obj.ContamParmIndices));
        end
        
        function parmList = RealsToParms(obj,Reals,~)
            % Return a list of the model parameter values retrieved from the reals used by fminsearch.
            % These come out in the order: TrueDists parms, prOutlier, Contam parms.
            parmList = zeros(obj.TotalNParms,1);
            for iCond = 1:obj.NConds
                theseIndices = obj.TrueParmIndices{iCond};
                parmList(theseIndices) = obj.TrueDists{iCond}.RealsToParms(Reals(theseIndices));
            end
            parmList(obj.prOutlierIndex) = NumTrans.Real2Bounded(obj.minPrOutlier,obj.maxPrOutlier,Reals(obj.prOutlierIndex));
            parmList(obj.ContamParmIndices) = obj.Contam.RealsToParms(Reals(obj.ContamParmIndices));
        end
        
    end % methods
    
end % clOutlierModel

