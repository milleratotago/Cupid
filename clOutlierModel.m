classdef clOutlierModel < handle
    % Model observed scores in NConds separate conditions as probability mixtures of scores
    % from a true distribution "contaminated" with some probability by an outlier process.
    %
    % Components of the model:
    %   TrueDist:   True distribution(s) of scores in the NConds conditions.
    %               This is either a single distribution, assumed to the the same across all conditions,
    %               or a cell array of NConds distributions that are allowed to differ across conditions.
    %   ContamDist: Distribution(s) of contamination values in the NConds conditions.
    %               This is either a single distribution, assumed to the the same across all conditions,
    %               or a cell array of NConds distributions that are allowed to differ across conditions.
    %   prOutlier:  Probability of an outlier.
    %               This is either a single value, assumed to be the same in all conditions,
    %               or a vector NConds values that are allowed to differ across conditions.
    %   ContamType: Indicates whether the contamination replaces the true score, adds to it, or multiplies it (see below)
    %               This is either a single value, assumed to be the same in all conditions,
    %               or a vector NConds values that are allowed to differ across conditions.
    %    
    % The parameters of an object of this class are assembled/stored in the order:
    %    TrueDist(s) parms, prOutlier(s), ContamDist(s) parms; (s) plurals refer to single/multiple conditions.
    
    properties(Constant)
        % ContamTypes are defined here.
        NTypes  = 3;  % N of types enumerated here.
        % Enum of different types of contamination:
        Shift   = 1;  % Contamination RV adds to true score
        Stretch = 2;  % Contamination RV multiplies times true score
        Replace = 3;  % Contamination RV replaces true score
        
        TypeLabel = {'Shift', 'Stretch', 'Replace'};  % Labels used for file naming, etc.
        
    end
    
    properties
        
        % "Input" properties defining the model:
        TrueDist    % Cell array of Cupid distribution(s) of true ("uncontaminated") scores in 1 or all NCond conditions.
        prOutlier   % Outlier probability, constant across conditions, or vector of probabilities.
        ContamDist  % Cell array of Cupid distribution(s) representing trial-to-trial variation in the contamination.
        ContamType  % Type of contamination: Shift, Stretch, or Replace, constant across conditions, or vector of those.
        
        % Optional input parameter pairs:
        %  'SearchOptions',optimset
        %  'SplinePDF',NSplinePDF
        %  'SplineCDF',NSplineCDF
        
        % Properties derived from the input properties:
        NConds     % N of conditions described by the model
        ObsDists   % Cell array of Cupid distributions of observed scores in the N conditions, including outliers.
        
        DefaultParmCodes
        SearchOptions
        SplinePDFs
        SplineCDFs

        % Counts
        nTrue        % length of TrueDist array
        nContam      % length of ContamDist array
        nprOutlier   % length of prOutlier vector
        nContamType  % length of ContamType vector
        
        % Indices of parameters within the overall parameter list:
        TrueParmIndices    % Cell array 1:nTrue of indices of the parameters for the true distribution(s)
        prOutlierIndices     % Vector 1:nprOutlier of index(es) of the prOutlier parameter(s)
        ContamParmIndices  % Cell array 1:nContam with indices of each outlier distribution's parameters
        TotalNParms

        % Vectors(1:NConds) of address pointers: Position iCond of each vector is the address
        % within the corresponding array of the position to be used for condition iCond.
        % Each vector is either 1,1,1,1... or 1:NConds depending on whether the referenced
        % variable is constant or varied across conditions.
        AddrTrue       %  Position within TrueDist to be used for condition iCond
        AddrContam     %  Position within ContamDist to be used for condition iCond
        AddrprOutlier  %  Position within prOutlier to be used for condition iCond
        AddrContamType       %  Position within ContamType to be used for condition iCond
        
    end
    
    properties (SetAccess = public)
        minPrOutlier, maxPrOutlier % Boundaries on prOutlier used in parameter searches.
    end
    
    methods
        
        function obj = clOutlierModel(TrueDist,prOutlier,ContamDist,ContamType,varargin)
            % Create object with specified values of input parameters.
            % TrueDist & ContamDist are either single Cupid distributions or cell arrays of Cupid distributions (1 per condition).
            % prOutlier & ContamType are either values or vectors of values (1 per condition).

            % Extract optional parameters:
            [obj.SearchOptions, varargin] = ExtractNameVali('SearchOptions',optimset,varargin);
            [obj.SplinePDFs, varargin] = ExtractNameVali({'SplinePDF','SplinePDFs'},0,varargin);
            [obj.SplineCDFs, varargin] = ExtractNameVali({'SplineCDF','SplineCDFs'},0,varargin);
            assert(numel(varargin)==0,'Unrecognized argument');

            % Make sure TrueDist and ContamDist are cell arrays (possibly just one element)
            if iscell(TrueDist)
                obj.TrueDist = TrueDist;
            else
                % Allow user to specify a single distribution simply as a Cupid
                % distribution rather than as a 1-position cell array.
                obj.TrueDist = {TrueDist};
            end
            if iscell(ContamDist)
                obj.ContamDist = ContamDist;
            else
                % Allow user to specify a single distribution simply as a Cupid
                % distribution rather than as a 1-position cell array.
                obj.ContamDist = {ContamDist};
            end

            % Count the number of inputs of each type:
            obj.nTrue = numel(TrueDist);
            obj.nContam = numel(ContamDist);
            obj.nprOutlier = numel(prOutlier);
            obj.nContamType = numel(ContamType);

            % Verify that all n's are 1 or max
            obj.NConds = max([obj.nTrue,obj.nContam,obj.nprOutlier,obj.nContamType]);
            NumsConsistent = ismember(obj.nTrue,[1,obj.NConds]) && ismember(obj.nContam,[1,obj.NConds]) && ...
                             ismember(obj.nprOutlier,[1,obj.NConds]) && ismember(obj.nContamType,[1,obj.NConds]);
            if ~NumsConsistent
                fprintf('nTrue = %d, nContam = %d, nprOutlier = %d, nContamType = %d\n',obj.nTrue,obj.nContam,obj.nprOutlier,obj.nContamType);
                error('Inconsistent numbers: all values other than 1 must be the same.');
            end

            % Verify appropriate values of prOutlier(s) & ContamType(s).
            for i=1:obj.nprOutlier
                assert(prOutlier(i) > 0 && prOutlier(i) < 1,'prOutlier(s) must be between 0 & 1');
            end
            for i=1:obj.nContamType
                assert(ContamType(i) >= 1 && ContamType(i) <= clOutlierModel.NTypes && ContamType(i) == floor(ContamType(i)),'Illegal ContamType(s)');
            end
            obj.prOutlier = prOutlier;
            obj.ContamType = ContamType;

            obj.minPrOutlier = min([0.001, prOutlier]);  % Defaults can be overridden
            obj.maxPrOutlier = max([0.999, prOutlier]);

            % Set up the vectors of address pointers.
            if obj.nTrue==1
                obj.AddrTrue = ones(obj.NConds,1);
            else
                obj.AddrTrue = 1:obj.NConds;
            end
            if obj.nContam==1
                obj.AddrContam = ones(obj.NConds,1);
            else
                obj.AddrContam = 1:obj.NConds;
            end
            if obj.nprOutlier==1
                obj.AddrprOutlier = ones(obj.NConds,1);
            else
                obj.AddrprOutlier = 1:obj.NConds;
            end
            if obj.nContamType==1
                obj.AddrContamType = ones(obj.NConds,1);
            else
                obj.AddrContamType = 1:obj.NConds;
            end

            obj.MakeObsDists;
            
            % Construct default parmcodes
            sParmCodes = '';
            for iTrue = 1:obj.nTrue
                sParmCodes = [sParmCodes obj.TrueDist{iTrue}.DefaultParmCodes]; %#ok<AGROW>
            end
            for iprOutlier=1:obj.nprOutlier
                sParmCodes = [sParmCodes 'r'];  %#ok<AGROW> prOutlier
            end
            for iContam = 1:obj.nContam
                sParmCodes = [sParmCodes obj.ContamDist{iContam}.DefaultParmCodes]; %#ok<AGROW>
            end
            obj.DefaultParmCodes = sParmCodes';  % ' to match with column vector whose size is checked in fminsearcharb
            
            % Compute the indices of all parameters within the overall parameter lists
            UsedParms = 0;
            obj.TrueParmIndices = cell(obj.nTrue,1);
            for iCond=1:obj.nTrue
                obj.TrueParmIndices{iCond} = (UsedParms+1):(UsedParms+obj.TrueDist{iCond}.NDistParms);
                UsedParms = UsedParms + obj.TrueDist{iCond}.NDistParms;
            end
%             obj.prOutlierIndices = zeros(obj.nprOutlier,1);
            obj.prOutlierIndices = (UsedParms+1):(UsedParms+obj.nprOutlier);
            UsedParms = UsedParms + obj.nprOutlier;
            obj.ContamParmIndices = cell(obj.nContam,1);
            for iCond=1:obj.nContam
                obj.ContamParmIndices{iCond} = (UsedParms+1):(UsedParms+obj.ContamDist{iCond}.NDistParms);
                UsedParms = UsedParms + obj.ContamDist{iCond}.NDistParms;
            end
            obj.TotalNParms = UsedParms;
            
        end
        
        function [] = MakeOneObsDist(obj,iCond)
            % Initialize the observed distribution for condition iCond the current parameter values.
            iTrue = obj.AddrTrue(iCond);
            iContam = obj.AddrContam(iCond);
            iprOutlier = obj.AddrprOutlier(iCond);
            iType = obj.AddrContamType(iCond);
            switch obj.ContamType(iType)
                case obj.Shift
                    obj.ObsDists{iCond} = ContamShift(obj.TrueDist{iTrue},obj.prOutlier(iprOutlier),obj.ContamDist{iContam});
                case obj.Stretch
                    obj.ObsDists{iCond} = ContamStretch(obj.TrueDist{iTrue},obj.prOutlier(iprOutlier),obj.ContamDist{iContam});
                case obj.Replace
                    obj.ObsDists{iCond} = Mixture(1-obj.prOutlier(iprOutlier),obj.TrueDist{iTrue},obj.ContamDist{iContam});
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
            Parm2RealFn = @obj.ParmsToReals;
            ErrFn = @MyErrFunc;
            StartingVals = ParmValues(obj);
            [EndingVals,fval,exitflag,output] = fminsearcharb(ErrFn,StartingVals,Real2ParmFn,Parm2RealFn,ParmCodes,obj.SearchOptions);
            for iCond=1:obj.NConds
                obj.ObsDists{iCond}.SkipImpossibleWarn = holdWarnStates(iCond);
            end
            obj.ResetParms(EndingVals);
            
            function thiserrval=MyErrFunc(realParms)
                ResetParms(obj,realParms)
                thiserrval = -LnLikelihood(obj,X);
            end
            
        end
        
        function ResetParms(obj,newparmvalues)
            % For each ObsDist distribution, find its adjusted parameters and then
            % call its reset parms function
            for iCond=1:obj.NConds
                iTrue = obj.AddrTrue(iCond);
                iContam = obj.AddrContam(iCond);
                iprOutlier = obj.AddrprOutlier(iCond);
                iType = obj.AddrContamType(iCond);
                newPrOutlier = newparmvalues(obj.prOutlierIndices(iprOutlier));
                obj.prOutlier(iprOutlier) = newPrOutlier;
                switch obj.ContamType(iType)
                    case obj.Replace
                        theseIndices = [obj.prOutlierIndices(iprOutlier), obj.TrueParmIndices{iTrue}, obj.ContamParmIndices{iContam}];
                        newparms4dist = newparmvalues(theseIndices);
                        newparms4dist(1) = 1 - newparms4dist(1);  % Replace uses mixture dist that wants PrInlier as 1st parm
                    case {obj.Shift, obj.Stretch}
                        theseIndices = [obj.TrueParmIndices{iTrue}, obj.prOutlierIndices(iprOutlier), obj.ContamParmIndices{iContam}];
                        newparms4dist = newparmvalues(theseIndices);
                    otherwise
                        error('Unrecognized ContamType.')
                end
                obj.ObsDists{iCond}.ResetParms(newparms4dist);
            end
        end
        
        function parms = ParmValues(obj)
            % Retrieve the current parameter values
            parms = zeros(1,obj.TotalNParms);
            for iCond=1:obj.nTrue
                theseIndices = obj.TrueParmIndices{iCond};
                parms(theseIndices) = obj.TrueDist{iCond}.ParmValues;
            end
            parms(obj.prOutlierIndices) = obj.prOutlier;
            for iCond=1:obj.nContam
                theseIndices = obj.ContamParmIndices{iCond};
                parms(theseIndices) = obj.ContamDist{iCond}.ParmValues;
            end
            
        end
        
        function realList = ParmsToReals(obj,Parms,~)
            % Return a list of the model parameter values transformed into reals as appropriate for fminsearch.
            % These come out in the order: TrueDist(s) parms, prOutlier(s), ContamDist(s) parms.
            realList = zeros(obj.TotalNParms,1);
            for iCond = 1:obj.nTrue
                theseIndices = obj.TrueParmIndices{iCond};
                realList(theseIndices) = obj.TrueDist{iCond}.ParmsToReals(Parms(theseIndices));
            end
            for iCond = 1:obj.nprOutlier
                realList(obj.prOutlierIndices(iCond)) = NumTrans.Bounded2Real(obj.minPrOutlier,obj.maxPrOutlier,Parms(obj.prOutlierIndices(iCond)));
            end
            for iCond = 1:obj.nContam
                theseIndices = obj.ContamParmIndices{iCond};
                realList(theseIndices) = obj.ContamDist{iCond}.ParmsToReals(Parms(theseIndices));
            end
        end
        
        function parmList = RealsToParms(obj,Reals,~)
            % Return a list of the model parameter values retrieved from the reals used by fminsearch.
            % These come out in the order: TrueDist(s) parms, prOutlier(s), ContamDist(s) parms.
            parmList = zeros(obj.TotalNParms,1);
            for iCond = 1:obj.nTrue
                theseIndices = obj.TrueParmIndices{iCond};
                parmList(theseIndices) = obj.TrueDist{iCond}.RealsToParms(Reals(theseIndices));
            end
            for iCond = 1:obj.nprOutlier
                parmList(obj.prOutlierIndices(iCond)) = NumTrans.Real2Bounded(obj.minPrOutlier,obj.maxPrOutlier,Reals(obj.prOutlierIndices(iCond)));
            end
            for iCond = 1:obj.nContam
                theseIndices = obj.ContamParmIndices{iCond};
                parmList(theseIndices) = obj.ContamDist{iCond}.RealsToParms(Reals(theseIndices));
            end
        end
        
    end % methods
    
end % clOutlierModel

%{
Old EstML variant that searched by iteration.
It did not work very well so I have not updated it to the latest version of clOutlierModel:

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
        

%}
