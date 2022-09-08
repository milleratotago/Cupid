classdef SpearKarDist < dContinuous
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
    
    methods (Static)
        
        % Function to produce monotonic probability estimates from observed
        % frequencies that may not yield monotonic probability estimates
        % due to binomial variability.
        % Provided by Rolf Ulrich & Karin Bausenhart, August 2021
        function [fiMono, nTrials, fi] = monotonize(n1, n2)
            % n1 and n2 are vectors of the same length as X.
            % At each k, n1(k) is the number of observations judged greater than X
            % and n2(k) is the number of observations judged less than X.
            
            nTrials = n1+n2;    % number of trials per comparison with X(k)
            fi = n2./nTrials;   % relative frequency of "less than X(k)" responses
            fiMono = fi;        % start value for monotonized data
            
            while any(fiMono(1:end-1) > fiMono(2:end))  % as long as there is any non-monotonicity
                i=1;
                while i <= length(fiMono) - 1      % check for all c-levels until the second-to-last
                    if fiMono(i) <= fiMono(i+1)    % if fi_c(i) <= fi_c(i+1)  (i.e., monotonous)
                        i=i+1;                              % do nothing but increase i
                    else                           % but if not monotonous
                        k=1;                                % start a counter k
                        while true
                            tempfi = sum(fiMono(i:i+k) .* nTrials(i:i+k)) ...  % compute temporary fi for fi_c(i) to fi_c(i+k)
                                ./ sum(nTrials(i:i+k));
                            if i+k+1 > length(fiMono); break                   % stop if i+k+1 < nclev
                            elseif fiMono(i+k+1) > tempfi; break               % stop when next fi would be larger than temporary mean (i.e., monotonous)
                            else
                                k=k+1;                                         % otherwise increase counter and start again, averaging over one more level
                            end
                        end
                        fiMono(i:i+k) = repmat(tempfi,1,k+1);                  % after exit, replace the fis from i:k with the monotonized fi
                        i=i+k+1;                                               % and just proceed at clevel i+k+1
                    end
                end
            end %
        end % monotonize fn
        
        function [Cs,monoCDFs] = Stretch01(Cs,monoCDFs)
            % Augment Cs and monoCDFs, if necessary, so that monoCDFs includes 0 and 1.
            % The choice of corresponding min and max Cs is somewhat arbitrary.
            if monoCDFs(1) > 0
                minC = Cs(1) - (Cs(2) - Cs(1));  % Make the overall minimum one "step" below the previous minimum
                Cs = [minC Cs];
                monoCDFs = [0 monoCDFs];
            end
            if monoCDFs(end) < 1
                maxC = Cs(end) + (Cs(end) - Cs(end-1));  % Make the overall maximum one "step" above the previous maximum
                Cs = [Cs maxC];
                monoCDFs = [monoCDFs 1];
            end
        end % Stretch01
        
    end % methods (Static)
    
    methods
        
        function obj=SpearKarDist(varargin)   % Constructor
            obj=obj@dContinuous('SpearKarDist');
            obj.ParmTypes = '';
            obj.DefaultParmCodes = '';
            obj.NDistParms = 0;
            % Handle constructor calls with different numbers of parameters.
            switch nargin
                case 0
                case 2
                    ResetParms(obj,varargin{1},varargin{2});
                otherwise
                    ME = MException('SpearKarDist:Constructor', ...
                        'Too many arguments passed to SpearKarDist constructor.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,Xs,CDFs)
            assert(issorted(Xs,'strictascend'),'Xs are not strictly increasing');
            assert(issorted(CDFs),'CDFs are not increasing');  % Not strict because CDF can stay constant as X increases.
            ClearBeforeResetParmsC(obj);
            obj.Xs = Xs;
            if numel(Xs) == numel(CDFs)
                obj.CDFs = CDFs;
            elseif numel(Xs) == numel(CDFs) + 2
                obj.CDFs = [0 CDFs 1];
            else
                error('Sizes of Xs and CDFs are not consistent.');
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
    
end  % class SpearKarDist

