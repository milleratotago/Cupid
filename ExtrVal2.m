classdef ExtrVal2 < dContinuous
    % Extreme Value Type II distribution with parameters mu, scale, shapek
    % Reference: JKN, Vol 2, page 3. (What they call epsilon I call mu, what they call theta I call scale)
    %  shapek is the shape parameter
    % Also mentioned by Luce (1986), page 508.
    
    properties(SetAccess = protected)
        mu, scale, shapek
    end
    
    methods
        
        function obj=ExtrVal2(varargin)
            obj=obj@dContinuous('ExtrVal2');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.SearchOptions.MaxFunEvals = 2000;
            obj.SearchOptions.MaxIter = 2000;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExtrVal2:Constructor', ...
                        'ExtrVal2 constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.scale = newparmvalues(2);
            obj.shapek = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newmu = ifelse(ParmCodes(1)=='f', obj.mu,obj.mu+obj.scale/100);
            newscale = ifelse(ParmCodes(2)=='f', obj.scale, 0.95*obj.scale);
            newshapek = ifelse(ParmCodes(3)=='f', obj.shapek, 1.05*obj.shapek);
            obj.ResetParms([newmu newscale newshapek]);
        end
        
        function []=ReInit(obj)
            assert(obj.scale>0,'ExtrVal2 scale must be > 0.');
            assert(obj.shapek>0,'ExtrVal2 shapek must be > 0.');
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            XX = (X(InBounds) - obj.mu) / obj.scale;
            XXP = XX.^(-obj.shapek);
            thispdf(InBounds) = obj.shapek*XXP.*exp(-XXP)./(XX*obj.scale);
            % for i=1:numel(X)
            %     if (X(i) > obj.LowerBound) && (X(i) < obj.UpperBound)
            %         XX = (X(i) - obj.mu) / obj.scale;
            %         XXP = XX^(-obj.shapek);
            %         thispdf(i) = obj.shapek*XXP*exp(-XXP)/(XX*obj.scale);
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            XX = (X(InBounds) - obj.mu) / obj.scale;
            XX = XX.^(-obj.shapek);
            thiscdf(InBounds) = exp(-XX);
            % for i=1:numel(X)
            %     if X(i) <= obj.LowerBound
            %     elseif X(i) >= obj.UpperBound
            %         thiscdf(i) = 1;
            %     else
            %         XX = (X(i) - obj.mu) / obj.scale;
            %         XX = XX^(-obj.shapek);
            %         thiscdf(i) = exp(-XX);
            %     end
            % end
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, ~, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            XX = -log(P);
            XX = XX.^(-1/obj.shapek);
            thisval = XX * obj.scale + obj.mu;
        end
        
    end  % methods
    
end  % class ExtrVal2



