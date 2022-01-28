classdef clOutlierModel2 < handle
    % Model observed scores in NConds separate conditions as probability mixtures of scores
    % from a true distribution "contaminated" with some probability by two outlier processes
    % (presumably unrealistically small and large values).
    %
    % LIMITATION: So far only "Replace" outliers are allowed.
    %
    % Components of the model:
    %   TrueDist:   True distribution(s) of scores in the NConds conditions.
    %               This is either a single distribution, assumed to the the same across all conditions,
    %               or a cell array of NConds distributions that are allowed to differ across conditions.
    %   ContamDist1: Distribution(s) of small contamination values in the NConds conditions.
    %                This is either a single distribution, assumed to the the same across all conditions,
    %                or a cell array of NConds distributions that are allowed to differ across conditions.
    %   prOutlier1:  Probability of a small outlier.
    %                This is either a single value, assumed to be the same in all conditions,
    %                or a vector NConds values that are allowed to differ across conditions.
    %   ContamDist2: Distribution(s) of large contamination values in the NConds conditions.
    %                This is either a single distribution, assumed to the the same across all conditions,
    %                or a cell array of NConds distributions that are allowed to differ across conditions.
    %   prOutlier2:  Probability of a large outlier.
    %                This is either a single value, assumed to be the same in all conditions,
    %                or a vector NConds values that are allowed to differ across conditions.
    %   ContamType: Indicates whether the contamination replaces the true score, adds to it, or multiplies it (see below)
    %                This is either a single value, assumed to be the same in all conditions,
    %                or a vector NConds values that are allowed to differ across conditions.
    %    
    % The parameters of an object of this class are assembled/stored in the order:
    %    TrueDist(s) parms, prOutlier1(s), prOutlier2(s), ContamDist1(s), ContamDist2(s) parms; (s) plurals refer to single/multiple conditions.
    
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
        prOutlier1  % Outlier probability, constant across conditions, or vector of probabilities.
        ContamDist1 % Cell array of Cupid distribution(s) representing trial-to-trial variation in the contamination.
        prOutlier2  % Outlier probability, constant across conditions, or vector of probabilities.
        ContamDist2 % Cell array of Cupid distribution(s) representing trial-to-trial variation in the contamination.
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
        nContam1     % length of ContamDist1 array
        nprOutlier1  % length of prOutlier1 vector
        nContam2     % length of ContamDist2 array
        nprOutlier2  % length of prOutlier2 vector
        nContamType  % length of ContamType vector
        
        % Indices of parameters within the overall parameter list:
        TrueParmIndices    % Cell array 1:nTrue of indices of the parameters for the true distribution(s)
        prOutlier1Indices  % Vector 1:nprOutlier1 of index(es) of the prOutlier1 parameter(s)
        prOutlier2Indices  % Vector 1:nprOutlier2 of index(es) of the prOutlier2 parameter(s)
        ContamParm1Indices % Cell array 1:nContam1 with indices of each outlier distribution's parameters
        ContamParm2Indices % Cell array 1:nContam2 with indices of each outlier distribution's parameters
        TotalNParms

        % Vectors(1:NConds) of address pointers: Position iCond of each vector is the address
        % within the corresponding array of the position to be used for condition iCond.
        % Each vector is either 1,1,1,1... or 1:NConds depending on whether the referenced
        % variable is constant or varied across conditions.
        AddrTrue        %  Position within TrueDist to be used for condition iCond
        AddrContam1     %  Position within ContamDist1 to be used for condition iCond
        AddrprOutlier1  %  Position within prOutlier1 to be used for condition iCond
        AddrContam2     %  Position within ContamDist2 to be used for condition iCond
        AddrprOutlier2  %  Position within prOutlier2 to be used for condition iCond
        AddrContamType  %  Position within ContamType to be used for condition iCond
        
    end
    
    properties (SetAccess = public)
        minprOutlier, maxprOutlier    % Boundaries on prOutlier used in parameter searches.
    end
    
    methods
        
        function obj = clOutlierModel2(TrueDist,prOutlier1,ContamDist1,prOutlier2,ContamDist2,ContamType,varargin)
            % Create object with specified values of input parameters.
            % TrueDist, ContamDist1, & ContamDist2 are either single Cupid distributions or cell arrays of Cupid distributions (1 per condition).
            % prOutlier1, prOutlier2 & ContamType are either values or vectors of values (1 per condition).

            % Extract optional parameters:
            [obj.SearchOptions, varargin] = ExtractNameVali('SearchOptions',optimset,varargin);
            [obj.SplinePDFs, varargin] = ExtractNameVali({'SplinePDF','SplinePDFs'},0,varargin);
            [obj.SplineCDFs, varargin] = ExtractNameVali({'SplineCDF','SplineCDFs'},0,varargin);
            assert(numel(varargin)==0,'Unrecognized argument');

            % Make sure TrueDist, ContamDist1, and ContamDist2 are cell arrays (possibly just one element)
            if iscell(TrueDist)
                obj.TrueDist = TrueDist;
            else
                % Allow user to specify a single distribution simply as a Cupid
                % distribution rather than as a 1-position cell array.
                obj.TrueDist = {TrueDist};
            end
            if iscell(ContamDist1)
                obj.ContamDist1 = ContamDist1;
            else
                % Allow user to specify a single distribution simply as a Cupid
                % distribution rather than as a 1-position cell array.
                obj.ContamDist1 = {ContamDist1};
            end
            if iscell(ContamDist2)
                obj.ContamDist2 = ContamDist2;
            else
                % Allow user to specify a single distribution simply as a Cupid
                % distribution rather than as a 1-position cell array.
                obj.ContamDist2 = {ContamDist2};
            end

            % Count the number of inputs of each type:
            obj.nTrue = numel(TrueDist);
            obj.nContam1 = numel(ContamDist1);
            obj.nprOutlier1 = numel(prOutlier1);
            obj.nContam2 = numel(ContamDist2);
            obj.nprOutlier2 = numel(prOutlier2);
            obj.nContamType = numel(ContamType);

            % Verify that all n's are 1 or max
            obj.NConds = max([obj.nTrue,obj.nContam1,obj.nprOutlier1,obj.nContam2,obj.nprOutlier2,obj.nContamType]);
            NumsConsistent = ismember(obj.nTrue,[1,obj.NConds]) && ismember(obj.nContamType,[1,obj.NConds]) && ...
                             ismember(obj.nprOutlier1,[1,obj.NConds]) && ismember(obj.nContam1,[1,obj.NConds]) && ...
                             ismember(obj.nprOutlier2,[1,obj.NConds]) && ismember(obj.nContam2,[1,obj.NConds]) ;
            if ~NumsConsistent
                fprintf('nTrue = %d, nContam1 = %d, nprOutlier1 = %d, nContam2 = %d, nprOutlier2 = %d, nContamType = %d\n',obj.nTrue,obj.nContam1,obj.nprOutlier1,obj.nContam2,obj.nprOutlier2,obj.nContamType);
                error('Inconsistent numbers: all values other than 1 must be the same.');
            end

            % Verify appropriate values of prOutlier1(s) & ContamType(s).
            for i=1:obj.nprOutlier1
                assert(prOutlier1(i) > 0 && prOutlier1(i) < 1,'prOutlier1(s) must be between 0 & 1');
            end
            for i=1:obj.nprOutlier2
                assert(prOutlier2(i) > 0 && prOutlier2(i) < 1,'prOutlier2(s) must be between 0 & 1');
            end
            for i=1:obj.nContamType
                assert(ContamType(i) >= 1 && ContamType(i) <= clOutlierModel2.NTypes && ContamType(i) == floor(ContamType(i)),'Illegal ContamType(s)');
            end
            obj.prOutlier1 = prOutlier1;
            obj.prOutlier2 = prOutlier2;
            obj.ContamType = ContamType;

            obj.minprOutlier = min([0.001, prOutlier1, prOutlier2]);  % Defaults can be overridden
            obj.maxprOutlier = max([0.999, prOutlier1, prOutlier2]);

            % Set up the vectors of address pointers.
            if obj.nTrue==1
                obj.AddrTrue = ones(obj.NConds,1);
            else
                obj.AddrTrue = 1:obj.NConds;
            end
            if obj.nContam1==1
                obj.AddrContam1 = ones(obj.NConds,1);
            else
                obj.AddrContam1 = 1:obj.NConds;
            end
            if obj.nContam2==1
                obj.AddrContam2 = ones(obj.NConds,1);
            else
                obj.AddrContam2 = 1:obj.NConds;
            end
            if obj.nprOutlier1==1
                obj.AddrprOutlier1 = ones(obj.NConds,1);
            else
                obj.AddrprOutlier1 = 1:obj.NConds;
            end
            if obj.nprOutlier2==1
                obj.AddrprOutlier2 = ones(obj.NConds,1);
            else
                obj.AddrprOutlier2 = 1:obj.NConds;
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
            for iprOutlier=1:obj.nprOutlier1
                sParmCodes = [sParmCodes 'r'];  %#ok<AGROW> prOutlier1
            end
            for iprOutlier=1:obj.nprOutlier2
                sParmCodes = [sParmCodes 'r'];  %#ok<AGROW> prOutlier2
            end
            for iContam = 1:obj.nContam1
                sParmCodes = [sParmCodes obj.ContamDist1{iContam}.DefaultParmCodes]; %#ok<AGROW>
            end
            for iContam = 1:obj.nContam2
                sParmCodes = [sParmCodes obj.ContamDist2{iContam}.DefaultParmCodes]; %#ok<AGROW>
            end
            obj.DefaultParmCodes = sParmCodes';  % ' to match with column vector whose size is checked in fminsearcharb
            
            % Compute the indices of all parameters within the overall parameter lists
            UsedParms = 0;
            obj.TrueParmIndices = cell(obj.nTrue,1);
            for iCond=1:obj.nTrue
                obj.TrueParmIndices{iCond} = (UsedParms+1):(UsedParms+obj.TrueDist{iCond}.NDistParms);
                UsedParms = UsedParms + obj.TrueDist{iCond}.NDistParms;
            end
            obj.prOutlier1Indices = (UsedParms+1):(UsedParms+obj.nprOutlier1);
            UsedParms = UsedParms + obj.nprOutlier1;
            obj.prOutlier2Indices = (UsedParms+1):(UsedParms+obj.nprOutlier2);
            UsedParms = UsedParms + obj.nprOutlier2;
            obj.ContamParm1Indices = cell(obj.nContam1,1);
            for iCond=1:obj.nContam1
                obj.ContamParm1Indices{iCond} = (UsedParms+1):(UsedParms+obj.ContamDist1{iCond}.NDistParms);
                UsedParms = UsedParms + obj.ContamDist1{iCond}.NDistParms;
            end
            obj.ContamParm2Indices = cell(obj.nContam2,1);
            for iCond=1:obj.nContam2
                obj.ContamParm2Indices{iCond} = (UsedParms+1):(UsedParms+obj.ContamDist2{iCond}.NDistParms);
                UsedParms = UsedParms + obj.ContamDist2{iCond}.NDistParms;
            end
            obj.TotalNParms = UsedParms;
            
        end
        
        function [] = MakeOneObsDist(obj,iCond)
            % Initialize the observed distribution for condition iCond the current parameter values.
            iTrue = obj.AddrTrue(iCond);
            iContam1 = obj.AddrContam1(iCond);
            iprOutlier1 = obj.AddrprOutlier1(iCond);
            iContam2 = obj.AddrContam2(iCond);
            iprOutlier2 = obj.AddrprOutlier2(iCond);
            iType = obj.AddrContamType(iCond);
            switch obj.ContamType(iType)
                case obj.Shift
                    error('Shift not supported yet.');
                case obj.Stretch
                    error('Shift not supported yet.');
                case obj.Replace
                    obj.ObsDists{iCond} = Mixture(1-obj.prOutlier1(iprOutlier1)-obj.prOutlier2(iprOutlier2),obj.TrueDist{iTrue},obj.prOutlier1(iprOutlier1),obj.ContamDist1{iContam1},obj.prOutlier2(iprOutlier2),obj.ContamDist2{iContam2});
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
            % X is a cell array, 1:NDataConds, of scores in the different conditions.
            % If the number of data conditions does not equal the number of model
            % conditions, then e.g. data conditions 1,2,3,4 correspond to model
            % conditions 1,2,1,2.
            X = EnsureCell(X);
            NDataConds = numel(X);
            perCondll = zeros(NDataConds,1);
            for iDataCond = 1:NDataConds
                iModelCond = JMod(iDataCond,obj.NConds);
                perCondll(iDataCond) = obj.ObsDists{iModelCond}.LnLikelihood(X{iDataCond});
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
                iContam1 = obj.AddrContam1(iCond);
                iprOutlier1 = obj.AddrprOutlier1(iCond);
                iContam2 = obj.AddrContam2(iCond);
                iprOutlier2 = obj.AddrprOutlier2(iCond);
                iType = obj.AddrContamType(iCond);
                newprOutlier1 = newparmvalues(obj.prOutlier1Indices(iprOutlier1));
                obj.prOutlier1(iprOutlier1) = newprOutlier1;
                newprOutlier2 = newparmvalues(obj.prOutlier2Indices(iprOutlier2));
                obj.prOutlier2(iprOutlier2) = newprOutlier2;
                switch obj.ContamType(iType)
                    case obj.Replace  %
                        % theseIndices picks out indices for TrueDist, prOutlier1, ContamParm1, ContamParm2
                        % Mixture distribution's ResetParms function wants all probabilities except the last, and PrInlier is added later
                        theseIndices = [obj.TrueParmIndices{iTrue}, obj.prOutlier1Indices(iprOutlier1), obj.ContamParm1Indices{iContam1}, obj.ContamParm2Indices{iContam2}];
                        newparms4dist = newparmvalues(theseIndices);
                        prOutlier1idx = obj.prOutlier1Indices(iprOutlier1);
                        tmpprOutlier1 = newparmvalues(prOutlier1idx);
                        prOutlier2idx = obj.prOutlier2Indices(iprOutlier2);
                        tmpprOutlier2 = newparmvalues(prOutlier2idx);
                        prInlier = 1 - tmpprOutlier1 - tmpprOutlier2;
                        newparms4dist = [prInlier, newparms4dist];  % Prepend the probability of the TrueDist
                    otherwise
                        error('Unrecognized/unsupported ContamType.')
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
            parms(obj.prOutlier1Indices) = obj.prOutlier1;
            parms(obj.prOutlier2Indices) = obj.prOutlier2;
            for iCond=1:obj.nContam1
                theseIndices = obj.ContamParm1Indices{iCond};
                parms(theseIndices) = obj.ContamDist1{iCond}.ParmValues;
            end
            for iCond=1:obj.nContam2
                theseIndices = obj.ContamParm2Indices{iCond};
                parms(theseIndices) = obj.ContamDist2{iCond}.ParmValues;
            end
            
        end
        
        function realList = ParmsToReals(obj,Parms,~)
            % Return a list of the model parameter values transformed into reals as appropriate for fminsearch.
            % These come out in the order: TrueDist(s) parms, prOutlier1(s), ContamDist1(s) parms.
            realList = zeros(obj.TotalNParms,1);
            for iCond = 1:obj.nTrue
                theseIndices = obj.TrueParmIndices{iCond};
                realList(theseIndices) = obj.TrueDist{iCond}.ParmsToReals(Parms(theseIndices));
            end
            for iCond = 1:obj.nprOutlier1
                realList(obj.prOutlier1Indices(iCond)) = NumTrans.Bounded2Real(obj.minprOutlier,obj.maxprOutlier,Parms(obj.prOutlier1Indices(iCond)));
            end
            for iCond = 1:obj.nprOutlier2
                realList(obj.prOutlier2Indices(iCond)) = NumTrans.Bounded2Real(obj.minprOutlier,obj.maxprOutlier,Parms(obj.prOutlier2Indices(iCond)));
            end
            for iCond = 1:obj.nContam1
                theseIndices = obj.ContamParm1Indices{iCond};
                realList(theseIndices) = obj.ContamDist1{iCond}.ParmsToReals(Parms(theseIndices));
            end
            for iCond = 1:obj.nContam2
                theseIndices = obj.ContamParm2Indices{iCond};
                realList(theseIndices) = obj.ContamDist2{iCond}.ParmsToReals(Parms(theseIndices));
            end
        end
        
        function parmList = RealsToParms(obj,Reals,~)
            % Return a list of the model parameter values retrieved from the reals used by fminsearch.
            % These come out in the order: TrueDist(s) parms, prOutlier1(s), ContamDist1(s) parms.
            parmList = zeros(obj.TotalNParms,1);
            for iCond = 1:obj.nTrue
                theseIndices = obj.TrueParmIndices{iCond};
                parmList(theseIndices) = obj.TrueDist{iCond}.RealsToParms(Reals(theseIndices));
            end
            for iCond = 1:obj.nprOutlier1
                parmList(obj.prOutlier1Indices(iCond)) = NumTrans.Real2Bounded(obj.minprOutlier,obj.maxprOutlier,Reals(obj.prOutlier1Indices(iCond)));
            end
            for iCond = 1:obj.nprOutlier2
                parmList(obj.prOutlier2Indices(iCond)) = NumTrans.Real2Bounded(obj.minprOutlier,obj.maxprOutlier,Reals(obj.prOutlier2Indices(iCond)));
            end
            for iCond = 1:obj.nContam1
                theseIndices = obj.ContamParm1Indices{iCond};
                parmList(theseIndices) = obj.ContamDist1{iCond}.RealsToParms(Reals(theseIndices));
            end
            for iCond = 1:obj.nContam2
                theseIndices = obj.ContamParm2Indices{iCond};
                parmList(theseIndices) = obj.ContamDist2{iCond}.RealsToParms(Reals(theseIndices));
            end
        end
        
    end % methods
    
end % clOutlierModel2

