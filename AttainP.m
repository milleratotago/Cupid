classdef AttainP < dContinuous % dEither  % Discrete not implemented yet. Maybe handle it by descending from dTransOf1 or dTransOf2 (MinBound)
    % AttainP(TrueDist,NullDist[,TailsArg]): Distribution of "p" values with true distribution tested against Null
    % TailsArg = -1 for lower tail p, 1 for upper tail p, or 2 for 2-tailed (default=2).');
    
    % This distribution is used to examine the distribution of attained p values
    % that will be obtained when data from a "True" data distribution
    % are judged for significance against a hypothesized "Null" distribution.
    %
    % For example, AttainP(Normal(1,1),Normal(0,1)) shows you the distribution of
    % p values that you would get if you took samples from a Normal(1,1) distribution
    % and tabulated the samples in terms of their attained 2-tailed p value's within the
    % Normal(0,1) distribution.
    %
    % Formally, for any p, 0<p<1, one establishes low and high cutoffs from the null
    % distribution such that Pr(X<=Low|Null) = p/2 and Pr(X>=Hi|Null) = p/2. The CDF
    % of p is then Pr(X<=Low|True) + Pr(X>=Hi|True).
    %
    % So far, the distribution only handles cases in which the True & Null distributions are both continuous.
    
    
    properties(SetAccess = protected)
        True, Null,
        % The next two booleans control options:
        TwoTailed,  % default=true,  so attained p value of True value is computed 2-tailed rather than 1-tailed.
        UpperTail   % default=false, relevant only when TwoTailed is false.
        % Then determines whether p value is computed from upper or lower tail of Null distribution.
    end
    
    properties(SetAccess = private)
        NTrueParms   % Holds the number of parameters in the True distribution.
        % These appear at the start of parameter lists.
    end
    
    methods
        
        function obj=AttainP(varargin)
            obj=obj@dContinuous('AttainP');
            switch nargin
                case 0
                case {2, 3}
                    Setup(obj,varargin(:));
                    ResetParms(obj,[obj.True.ParmValues obj.Null.ParmValues]);
                otherwise
                    ME = MException('AttainP:Constructor', ...
                        'AttainP constructor needs 0, 2, or 3 arguments.');
                    throw(ME);
            end
        end
        
        function Setup(obj,s)
            obj.True = s{1};
            obj.Null = s{2};
            obj.DistType = obj.True.DistType;
            parmspassed = numel(s);
            if (parmspassed == 2) || (s{3}==2)
                obj.TwoTailed = true;
                obj.UpperTail = false;
            else
                switch s{3}
                    case -1
                        obj.TwoTailed = false;
                        obj.UpperTail = false;
                    case 1
                        obj.TwoTailed = false;
                        obj.UpperTail = true;
                    otherwise
                        ME = MException('AttainP:Setup', ...
                            'AttainP NTails argument must be -1, 1, or 2 (default=2).');
                        throw(ME);
                end
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = ['AttainP(' obj.True.StringName ',' obj.Null.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.Initialized = false;
            obj.True.ResetParms(newparmvalues(1:obj.True.NDistParms));
            obj.Null.ResetParms(newparmvalues(obj.True.NDistParms+1:end));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.True.PerturbParms(ParmCodes(1:obj.True.NDistParms));
            obj.Null.PerturbParms(ParmCodes(obj.True.NDistParms+1:end));
            obj.ResetParms([obj.True.ParmValues obj.Null.ParmValues]);
        end
        
        function []=ReInit(obj)
            obj.NTrueParms = obj.True.NDistParms;
            obj.NDistParms = obj.Null.NDistParms + obj.True.NDistParms;
            s = obj.Null.DefaultParmCodes;
            s(1:end) = 'f';
            obj.DefaultParmCodes = [obj.True.DefaultParmCodes s];
            obj.Initialized = true;
            obj.DistType = obj.True.DistType;
            obj.NValues = obj.True.NValues;
            %             if (obj.True.DistType == 'c') && (obj.Null.DistType == 'c')
            %                 obj.DistType = 'c';
            if obj.TwoTailed
                LB1 = obj.Null.CDF(obj.True.LowerBound);
                LB2 = 1 - obj.Null.CDF(obj.True.UpperBound);
                if LB1 <= LB2
                    obj.LowerBound = LB1;
                else
                    obj.LowerBound = LB2;
                end
                NullMed = obj.Null.Median;
                if (obj.True.LowerBound <= NullMed) && (obj.True.UpperBound >= NullMed)
                    obj.UpperBound = 1;
                elseif obj.True.LowerBound > NullMed
                    obj.UpperBound = 2*(1-obj.Null.CDF(obj.True.LowerBound));
                    if obj.True.UpperBound < NullMed
                        obj.UpperBound = 2*obj.Null.CDF(obj.True.UpperBound);
                    end
                end  % TwoTailed
            else   % Not TwoTailed
                if obj.UpperTail
                    obj.LowerBound = 1 - obj.Null.CDF(obj.True.UpperBound);
                    obj.UpperBound = 1 - obj.Null.CDF(obj.True.LowerBound);
                else
                    obj.LowerBound = obj.Null.CDF(obj.True.LowerBound);
                    obj.UpperBound = obj.Null.CDF(obj.True.UpperBound);
                end
            end  % Not TwoTailed
            %                 % end of "if 'c'"
            %             else
            %                 msgID = 'AttainP:IllegalDist';
            %                 msg = 'Unable to handle this class of distributions.';
            %                 ThisException = MException(msgID,msg);
            %                 throw(ThisException);
            %             end
            if obj.NameBuilding
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            % From Ulrich Miller (2016) p-curve paper
            for iel=1:numel(X)
                if InBounds(iel)
                    if obj.TwoTailed
                        ThisXLo = obj.Null.InverseCDF(X(iel)/2);
                        ThisXHi = obj.Null.InverseCDF(1-X(iel)/2);
                        RatioLo = obj.True.PDF(ThisXLo) / obj.Null.PDF(ThisXLo);
                        RatioHi = obj.True.PDF(ThisXHi) / obj.Null.PDF(ThisXHi);
                        thispdf(iel) = 0.5*(RatioLo + RatioHi);
                    elseif obj.UpperTail
                        ThisX = obj.Null.InverseCDF(1-X(iel));
                        thispdf(iel) = obj.True.PDF(ThisX) / obj.Null.PDF(ThisX);
                    else
                        ThisX = obj.Null.InverseCDF(X(iel));
                        thispdf(iel) = obj.True.PDF(ThisX) / obj.Null.PDF(ThisX);
                    end
                end
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    if obj.TwoTailed
                        thiscdf(iel) = obj.True.CDF(obj.Null.InverseCDF(X(iel)/2)) + 1 - obj.True.CDF(obj.Null.InverseCDF(1-X(iel)/2));
                    elseif obj.UpperTail
                        thiscdf(iel) = 1 - obj.True.CDF(obj.Null.InverseCDF(1-X(iel)));
                    else
                        thiscdf(iel) = obj.True.CDF(obj.Null.InverseCDF(X(iel)));
                    end
                end
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.True.ParmsToReals(Parms(1:obj.True.NDistParms)) obj.Null.ParmsToReals(Parms(obj.True.NDistParms+1:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.True.RealsToParms(Reals(1:obj.True.NDistParms)) obj.Null.RealsToParms(Reals(obj.True.NDistParms+1:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.True.ParmValues obj.Null.ParmValues];
        end
        
        % function thisval = nIthValue(iVal)
        %     thisval = 1;  % Dummy
        % end
        
        % function thisval = NearestLegal(X)
        %     thisval = 1;  % Dummy
        % end
        
        % function thisval = LegalValue(X)
        %     thisval = false;  % Dummy
        % end
        
        function x = XsToPlot(obj)
            x = XsToPlot(obj.True);
            for i=1:numel(x)
                if obj.TwoTailed
                    if x(i) < 0
                        x(i) = CDF(obj.Null,x(i)) + 1 - CDF(obj.Null,-x(i));
                    else
                        x(i) = CDF(obj.Null,-x(i)) + 1 - CDF(obj.Null,x(i));
                    end
                elseif obj.UpperTail
                    x(i) = 1 - CDF(obj.Null,x(i));
                else
                    x(i) = CDF(obj.Null,x(i));
                end
            end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            r = Random(obj.True,varargin{:});
            thisval=zeros(varargin{:});
            for i=1:numel(thisval)
                if obj.TwoTailed
                    if r(i) < 0
                        thisval(i) = CDF(obj.Null,r(i)) + 1 - CDF(obj.Null,-r(i));
                    else
                        thisval(i) = CDF(obj.Null,-r(i)) + 1 - CDF(obj.Null,r(i));
                    end
                elseif obj.UpperTail
                    thisval(i) = 1 - CDF(obj.Null,r(i));
                else
                    thisval(i) = CDF(obj.Null,r(i));
                end
            end
        end
        
    end  % methods
    
end  % class AttainP


