classdef dMATLABc < dContinuous
    % dMATLABc(pd)   pd is any continuous probability distribution object in MATLAB
    % This class serves as an interface between Cupid and MATLAB's generic distributions.
    % For examples, see Demo_dMATLAC.m

    % To use some Cupid functionality (e.g., parameter estimation), the user must provide
    % the following 3 pieces of information about each parameter of the distribution:
    %    integer or real?
    %    minimum legal value (may be -inf)
    %    maximum legal value (may be +inf)
    % This information is provided via three vectors, each with length equal to the number
    %   of distribution parameters:
    %   RealInt: character array of 'i' or 'r' values.
    %   MinVal and MaxVal: real arrays of the min & max values, which may be -inf and +inf
    %    for parameters unbounded in either or both directions.
    %   Note: Parameters bounded only from above are not supported.
    % The three vectors can be provided either:
    %    when the constructor is called to create the object, or
    %    after construction, by calling SetParmInfo(obj,RealInt,MinVal,MaxVal).
    
    properties(SetAccess = protected)
        pd    % The user is not allowed to alter pd after the object is created.
        ParmBoundClass  % Array 1..NDistParms of 0/1/2, where
                        % 0 = unbounded; no need for NumTrans
                        % 1 = bounded below; use GT2Real and Real2GT
                        % 3 = bounded below & above; use Bounded2Real and Real2Bounded
        MinParmVal, MaxParmVal  % Arrays of minimum and maximum legal parameter values
    end
    
    methods
        
        function obj=dMATLABc(passpd,varargin)
            obj=obj@dContinuous('dMATLABc');
            obj.pd = passpd;
            % Next line removes hyphens & spaces from e.g., Birnbaum-Saunders DistributionName
            obj.FamilyName = regexprep(passpd.DistributionName,'[- ]','');
            if isprop(obj.pd,'NumParameters')
                obj.NDistParms = obj.pd.NumParameters;
            else
                obj.NDistParms = 0;  % e.g., kernel density distributions have no parameters.
            end
            obj.DefaultParmCodes = repmat('r',1,obj.NDistParms);

            switch numel(varargin)
                case 0
                case 3
                    obj.SetParmInfo(varargin{:});
                otherwise
                    ME = MException('dMATLABc:Constructor','dMATLABc constructor needs 1 or 4 arguments.');
                    throw(ME);
            end
            obj.ReInit;
        end
        
        function SetParmInfo(obj,RealInt,MinVal,MaxVal)
            obj.ParmTypes = RealInt;
            obj.DefaultParmCodes = RealInt;
            obj.MinParmVal = MinVal;
            obj.MaxParmVal = MaxVal;
            obj.ParmBoundClass = zeros(obj.NDistParms,1);
            for iParm=1:obj.NDistParms
                if ~isinf(obj.MinParmVal(iParm))
                   obj.ParmBoundClass(iParm) = 1;
                end
                if ~isinf(obj.MaxParmVal(iParm))
                   obj.ParmBoundClass(iParm) = obj.ParmBoundClass(iParm) + 2;
                end
                assert(~(obj.ParmBoundClass(iParm)==2),'Bounding of parameters only from above is not supported.');
            end
        end

        function parmvals = ParmValues(obj)
            parmvals = zeros(1,obj.NDistParms);
            for i=1:obj.NDistParms
                try  % PiecewiseLinear fails because the parameters are the vectors x and Fx.
                    parmvals(i) = obj.pd.(obj.pd.ParameterNames{i});
                catch
                    parmvals(i) = nan;
                end
            end
        end
        
        function BuildMyName(obj)
            parmvals = obj.ParmValues;
            % parmvalstrs = num2str(parmvals);
            parmvalstrs = regexprep(num2str(parmvals),'\s+',',');  % Comma separated numbers
            obj.StringName = [obj.FamilyName '(' parmvalstrs ')'];
        end

        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            parmnamesandvals = cell(1,2*obj.NDistParms);
            for i=1:obj.NDistParms
                parmnamesandvals{2*i-1} = obj.pd.ParameterNames{i};
                parmnamesandvals{2*i} = newparmvalues(i);
            end
            obj.pd = makedist(obj.FamilyName,parmnamesandvals{:});
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newparms = obj.ParmValues;
            for iParm=1:obj.NDistParms
                if ~(ParmCodes(iParm)=='f')
                    switch obj.ParmBoundClass(iParm)
                        case {0, 1}
                            newparms(iParm) = newparms(iParm) * 1.01 + 0.1;
                        case 3
                            if newparms(iParm) <= (obj.MinParmVal(iParm)+obj.MaxParmVal(iParm)) / 2
                                newparms(iParm) = newparms(iParm) + 0.1*(obj.MaxParmVal(iParm) - newparms(iParm));
                            else
                                newparms(iParm) = newparms(iParm) - 0.1*(newparms(iParm) - obj.MinParmVal(iParm));
                            end
                    end  % switch
                end  % if
            end  % for
            obj.ResetParms(newparms);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
            obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = Parms;
            for iParm=1:obj.NDistParms
                switch obj.ParmBoundClass(iParm)
                    case 1
                        Reals(iParm) = NumTrans.GT2Real(obj.MinParmVal(iParm),Parms(iParm));
                    case 3
                        Reals(iParm) = NumTrans.Bounded2Real(obj.MinParmVal(iParm),obj.MaxParmVal(iParm),Parms(iParm));
                end
            end
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = Reals;
            for iParm=1:obj.NDistParms
                switch obj.ParmBoundClass(iParm)
                    case 1
                        Parms(iParm) = NumTrans.Real2GT(obj.MinParmVal(iParm),Reals(iParm));
                    case 3
                        Parms(iParm) = NumTrans.Real2Bounded(obj.MinParmVal(iParm),obj.MaxParmVal(iParm),Reals(iParm));
                end
            end
        end
        
        function thispdf=PDF(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.pd.pdf(X(InBounds));
        end
        
        function thiscdf=CDF(obj,X)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = obj.pd.cdf(X(InBounds));
        end
        
        function thisval=InverseCDF(obj,P)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.pd.icdf(P);
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.pd.random(varargin{:});
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.pd.mean;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.pd.var;
        end
        
        function thisval=SD(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.pd.std;
        end
        
%         function s=EstML2(obj,x)  % Use of MATLAB fitdist after object creation is not supported
%             newpd = fitdist(x,obj.FamilyName)  % There are many more parameters for this.
%             obj.pd = newpd;
% This is not quite right because the parameters are not reset (e.g., bounds).
% It is better anyway to let the user use fitdist as desired and make a new dMATLABc
% from the result.
%             BuildMyName(obj);
%             s=obj.StringName;
%         end
        
    end  % methods
    
end  % class dMATLABc

