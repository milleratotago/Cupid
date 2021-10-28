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

    methods
        
        function obj = clOutlierModel(TrueDists,prOutlier,Contam,ConType,varargin)
            % Create object with specified values of input parameters.
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
        
        function [] = MakeObsDists(obj)
            % Initialize the observed distributions with the current parameter values.
            obj.ObsDists = cell(obj.NConds,1);
            for iCond=1:obj.NConds
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
            Real2ParmFn = @obj.RealsToParms;
            Parm2ReamFn = @obj.ParmsToReals;
            ErrFn = @MyErrFunc;
            StartingVals = ParmValues(obj);
            [EndingVals,fval,exitflag,output] = fminsearcharb(ErrFn,StartingVals,Real2ParmFn,Parm2ReamFn,ParmCodes,obj.SearchOptions);
            obj.ResetParms(EndingVals);
           
             function thiserrval=MyErrFunc(realParms)
                ResetParms(obj,realParms)
                thiserrval = -LnLikelihood(obj,X);
            end
            
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
            realList(obj.prOutlierIndex) = NumTrans.Bounded2Real(0,1,Parms(obj.prOutlierIndex));
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
            parmList(obj.prOutlierIndex) = NumTrans.Real2Bounded(0,1,Reals(obj.prOutlierIndex));
            parmList(obj.ContamParmIndices) = obj.Contam.RealsToParms(Reals(obj.ContamParmIndices));
        end

    end % methods
    
end % clOutlierModel

