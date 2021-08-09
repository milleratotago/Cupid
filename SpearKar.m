classdef SpearKar < dContinuous
    % Spearman-Kaerber distribution
    % This distribution is characterized by a discrete set of k X's: X(1)...X(k) and
    %   k associated CDF values such that Pr(x<X(i) = CDF(i).
    % It is assumed that the CDF increases linearly from X(i) to X(i+1)
    % Based on that assumption, PDFs are computed for all points between
    %  the minimum X(1) and the maximum X(k)

    properties
        Xs   % values covering the range of possible X's.  X(1) is min and X(end) is max
        CDFs % CDFs at each X, which must start at 0 and increase to 1
        PDFs % from Xi to X_{i+1}  (there are only k-1 of these)
    end

    methods
        
        function obj=SpearKar(varargin)   % Constructor
            obj=obj@dContinuous('SpearKar');
            obj.ParmTypes = '';
            obj.DefaultParmCodes = '';
            obj.NDistParms = 0;
            % Handle constructor calls with different numbers of parameters.
            switch nargin
                case 0
                case 2
                    ResetParms(obj,varargin{1},varargin{2});
                otherwise
                    ME = MException('SpearKar:Constructor', ...
                        'Too many arguments passed to SpearKar constructor.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,Xs,CDFs)
            assert(issorted(Xs,'strictascend'),'Xs are not strictly increasing');
            assert(issorted(CDFs,'strictascend'),'CDFs are not strictly increasing');
            ClearBeforeResetParmsC(obj);
            obj.Xs = Xs;
            if numel(Xs) == numel(CDFs)
                obj.CDFs = CDFs;
            elseif numel(Xs) == numel(CDFs) + 2
                obj.CDFs = [0 CDFs 1];
            else
                error('Sizes of Xs and CDFs do not match');
            end
            assert(obj.CDFs(1)==0,'first CDF is not zero');
            assert(obj.CDFs(end)==1,'last CDF is not one');
            ReInit(obj);
        end
        
%       function PerturbParms(obj,ParmCodes)
%           if ~(ParmCodes(1) == 'f')
%               obj.ResetParms(obj.df+1);
%           end
%       end
        
        function Reals = ParmsToReals(obj,Parms,~)  %#ok<INUSL>
            Reals = Parms;  % Required but unused method
        end
        
        function Parms = RealsToParms(obj,Reals,~) %#ok<INUSL>
            Parms = Reals;  % Required but unused method
        end

        function []=ReInit(obj)
            obj.LowerBound = obj.Xs(1);
            obj.UpperBound = obj.Xs(end);
            obj.PDFs = obj.Xs(1:end-1);  % just making a vector of the right length
            for i=1:numel(obj.PDFs)
                range = obj.Xs(i+1) - obj.Xs(i);
                prob  = obj.CDFs(i+1) - obj.CDFs(i);
                obj.PDFs(i) = prob / range;
            end
            obj.Initialized = true;  % Needed so that CDF can be called
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            thispdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            for i=1:numel(X)
               if InBounds(i)
                  iabove = find(X(i)<=obj.Xs,1);  % find first obj.Xs < X(i)
                  iabove = max(iabove,2);  % prevent error when X == obj.Xs(1)
                  thispdf(i) = obj.PDFs(iabove-1);
               end
            end
        end
        
        function thiscdf=CDF(obj,X)
            thiscdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(X>=obj.UpperBound) = 1;
            for i=1:numel(X)
               if InBounds(i)
                  iabove = find(X(i)<=obj.Xs,1);  % find first obj.Xs < X(i)
                  iabove = max(iabove,2);  % prevent error when X == obj.Xs(1)
                  dist = X(i) - obj.Xs(iabove-1);
                  pdfhere = obj.PDFs(iabove-1);
                  thiscdf(i) = obj.CDFs(iabove-1) + dist*pdfhere;
               end
            end
        end
        
        function thisval=Mean(obj)
            thisval = obj.RawMoment(1);
        end
        
        function thisval=Variance(obj)
            thisval = obj.RawMoment(2) - obj.RawMoment(1)^2;
        end
        
        function thisval=RawMoment(obj,r)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            total = 0;
            for i=2:numel(obj.Xs)
                total = total + (obj.CDFs(i) - obj.CDFs(i-1)) * (obj.Xs(i)^(r+1) - obj.Xs(i-1)^(r+1)) / (obj.Xs(i) - obj.Xs(i-1));
            end
            thisval = total / (r+1);
        end
        
    end  % methods
    
end  % class SpearKar

