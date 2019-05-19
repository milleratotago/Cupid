classdef F < dContinuous
    % F(dfnum,dfdenom)

    % WARNINGS:
    %   Changes to this distribution also affect some other distributions using it.
    %   There are numerical errors in computing moments, especially bad with small dfs.
    
    properties(SetAccess = protected)
        dfnum, dfdenom,
        LnDF1, LnDF2, LnGamHalfSum, LnGamHalf1, LnGamHalf2,
        ComputeBounds, BetaGenerator
    end
    
    methods
        
        function obj=F(varargin)
            obj=obj@dContinuous('F');
            obj.ParmTypes = 'ii';
            obj.DefaultParmCodes = 'ii';
            obj.NDistParms = 2;
            obj.BetaGenerator = Beta;
            obj.ComputeBounds = true;  % Estimate bounds to achieve CDF's of CDFNearlyZero and CDFNearlyOne
            obj.CDFNearlyOne = 1 - 1e-10;  % Try to get very close to end of long positive tail.
            %            obj.IntegrateOverP = true;  % This did not help with moment computations for F(1,6)
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('F:Constructor', ...
                        'F constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.dfnum = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.dfdenom = VerifyIntegerGE(obj,1,newparmvalues(2));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            Newnum = ifelse(ParmCodes(1)=='f', obj.dfnum, obj.dfnum+1);
            Newdenom = ifelse(ParmCodes(2)=='f', obj.dfdenom, obj.dfdenom+1);
            obj.ResetParms([Newnum Newdenom]);
        end
        
        function []=ReInit(obj)
            obj.LnDF1 = log(obj.dfnum);
            obj.LnDF2 = log(obj.dfdenom);
            obj.LnGamHalfSum = gammaln( (obj.dfnum+obj.dfdenom)/2 );
            obj.LnGamHalf1 = gammaln( obj.dfnum / 2 );
            obj.LnGamHalf2 = gammaln( obj.dfdenom / 2 );
            ResetParms(obj.BetaGenerator,[obj.dfnum/2 obj.dfdenom/2]);
            obj.Initialized = true;
            obj.LowerBound = 1e-15;
            obj.UpperBound = 1000;
            if obj.ComputeBounds
                obj.LowerBound = 0; % 1e-10;
                obj.UpperBound = 10000;
                %                 obj.UpperBound = 1.5;
                %                 while CDF(obj,obj.UpperBound) < obj.CDFNearlyOne
                %                     obj.UpperBound = 2*obj.UpperBound;
                %                 end
                obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
                obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.GT2Real(1,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2GT(1,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            % Based on F-distribution in CRC handbook of tables:
            Lnthispdf = obj.LnGamHalfSum - obj.LnGamHalf1 - obj.LnGamHalf2 ...
                + obj.LnDF1 * (obj.dfnum / 2) ...
                + obj.LnDF2 * (obj.dfdenom / 2) ...
                + log(X(InBounds)) * (obj.dfnum - 2) / 2 ...
                - log(obj.dfdenom + obj.dfnum * X(InBounds)) * (obj.dfnum + obj.dfdenom) / 2;
            thispdf(InBounds) = exp(Lnthispdf);
            % for i=1:numel(X)
            %     if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %         % Based on F-distribution in CRC handbook of tables:
            %         Lnthispdf = obj.LnGamHalfSum - obj.LnGamHalf1 - obj.LnGamHalf2 ...
            %             + obj.LnDF1 * (obj.dfnum / 2) ...
            %             + obj.LnDF2 * (obj.dfdenom / 2) ...
            %             + log(X(i)) * (obj.dfnum - 2) / 2 ...
            %             - log(obj.dfdenom + obj.dfnum * X(i)) * (obj.dfnum + obj.dfdenom) / 2;
            %         thispdf(i) = exp(Lnthispdf);
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = 1 - betainc(obj.dfdenom ./ (obj.dfdenom+obj.dfnum*X(InBounds)),obj.dfdenom/2,obj.dfnum/2);
            % for i=1:numel(X)
            %     if InBounds(i)
            %         thiscdf(i) = 1 - betainc(obj.dfdenom/(obj.dfdenom+obj.dfnum*X(i)),obj.dfdenom/2,obj.dfnum/2);
            %     end
            % end
        end
        
        function thisval=Mean(obj) % JKB, Vol 2, p 326
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.dfdenom > 2
                thisval = obj.dfdenom / (obj.dfdenom - 2);
            else
                thisval = inf;
            end
        end
        
        function thisval=InverseCDF(obj,P)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            BetaX = InverseCDF(obj.BetaGenerator,P);
            thisval = obj.dfdenom / obj.dfnum * BetaX ./ (1 - BetaX);
        end
        
        function thisval=Variance(obj)   % Wikipedia & JKB, Vol 2, p 326
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.dfdenom <= 4
                thisval = inf;  %  'Requested infinite variance of F distribution.'
            else
                thisval = 2 * obj.dfdenom^2 * (obj.dfnum + obj.dfdenom - 2) / (  obj.dfnum * (obj.dfdenom - 2)^2 * (obj.dfdenom - 4)  );
            end
        end
        
        function thisval=Skewness(obj)   % Wikipedia
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.dfdenom <= 6
                thisval = inf;  %  'Requested infinite skewness of F distribution.'
            else
                thisval = (2*obj.dfnum+obj.dfdenom-2)*sqrt(8*(obj.dfdenom-4)) / ...
                    ( (obj.dfdenom-6)*sqrt(obj.dfnum*(obj.dfnum+obj.dfdenom-2)) );
            end
        end
        
        %        function thisval=Kurtosis(obj)   % Dec 2016 Wikipedia gives this function for excess kurtosis, but I think it is wrong.
%             if ~obj.Initialized
%                 error(UninitializedError(obj));
%             end
        %            if obj.dfdenom <= 8
        %                thisval = inf;  %  'Requested infinite kurtosis of F distribution.'
        %            else
        %                thisval = 3 + 12*obj.dfnum*(5*obj.dfdenom-22)*(obj.dfnum+obj.dfdenom-2) + (obj.dfdenom-4)*(obj.dfdenom-2)^2 / ...
        %                          ( obj.dfnum*(obj.dfdenom-6)*(obj.dfdenom-8)*(obj.dfnum+obj.dfdenom-2) );
        %            end
        %        end
        
        function thisval=RawMoment(obj,I)   % JKB, Vol 2, p 325
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if I == 0
                thisval = 1;
            elseif I == 1
                thisval = Mean(obj);
            elseif I >= obj.dfdenom/2
                thisval = inf;
                warning(['Requested infinite raw moment ' num2str(I) ' of F distribution with dfdenom = ' num2str(obj.dfdenom)])
            else
                %                dfnumProd = obj.dfnum;
                %                for J = 1:I-1
                %                    dfnumProd = dfnumProd * (obj.dfnum + 2 * J);
                %                end
                %                dfdenomProd = 1;
                %                for J = 1:I
                %                    dfdenomProd = dfdenomProd * (obj.dfdenom - 2 * J);
                %                end
                %                thisval = (obj.dfdenom/obj.dfnum)^I * dfnumProd / dfdenomProd;
                thisval = (obj.dfdenom/obj.dfnum)^I * gamma(obj.dfnum/2 + I) / gamma(obj.dfnum/2) * gamma(obj.dfdenom/2 - I) / gamma(obj.dfdenom/2);
            end
        end
        
        function thisval=Random(obj,varargin)  % From Beta: Devroye p 430
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            X = Random(obj.BetaGenerator,varargin{:});
            thisval = obj.dfdenom / obj.dfnum * X ./ (1 - X);
        end
        
        function []=EstimatesFromMoments(obj,PassMoments, ParmCodes)
            HoldComputeBounds = obj.ComputeBounds;
            obj.ComputeBounds = numel(PassMoments) > 2;
            EstimatesFromMoments@dContinuous(PassMoments, ParmCodes);
            obj.ComputeBounds = HoldComputeBounds;
        end
        
        function thisval=MGF(obj,~)
            warning('MGF of F distribution does not exist.');
            thisval = nan;
        end
        
    end  % methods
    
end  % class F





