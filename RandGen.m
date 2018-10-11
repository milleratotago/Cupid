classdef RandGen < handle
    % Tools for multivariate random number generation.
    
% Copyright (C) 2018 Jeffrey Owen Miller
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program (00License.txt). If not, see 
%     <http://www.gnu.org/licenses/>.
%
    
    properties(SetAccess = public)  % These properties can be changed without restriction
        WantHistograms    % Set to true to request plotting histograms for every random sample
        WantScattergrams  % Set to false to request plotting scattergrams for every random sample
    end
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        NRVs
        RVs    % RVs is a cell array of random variables
               % RVs{1}, RVs{2}, etc are individual random variables, each one of which is descended from dContinuous.
        NSteps, Pctiles, ZPctiles, RVPctiles  % Used when approximating distributions to adjust RhoControllers
        RhoControllers, TargetCorrs, RhoMax, RhoMin
    end
    
    methods (Static)
        
        function rands = GenRands(RVs,RhoControllers,NCases,varargin)
            % This function generates correlated random variables, and it may be
            % the only function of this class that you will ever need to call.
            % The disadvantage of this function is that it saves no information,
            % so everything has to be recomputed every time you call it.
            % If you want to generate random numbers from the same distribution
            % repeatedly, you should make your own an object.
            obj = RandGen(RVs,RhoControllers,varargin{:});
            rands = obj.Rands(NCases);
        end
        
    end % Static methods
    
    methods
        
        function obj=RandGen(varargin)   % Constructor
            % The constructor may be called with no arguments or with arguments
            % in the ordered needed for Init, which is (RVs,RhoControllers,varargin).
            obj.WantHistograms = false;
            obj.WantScattergrams = false;
            if numel(varargin)>=2
                Init(obj,varargin{1},varargin{2},varargin{3:end});
            end
        end % Constructor
        
        function Init(obj,RVs,RhoControllers,varargin)
            obj.RVs = RVs;
            obj.NRVs = numel(RVs);
            obj.RhoControllers = RhoControllers;
            try
                [obj.WantHistograms, varargin] = ExtractNamei({'Histogram','Histograms'},varargin);
            catch ME
                disp('Error handling varargin with ExtractName');
                disp('For RandGen, you also need ExtractNameVal from https://github.com/milleratotago/ExtractNameVal ');
                rethrow(ME);
            end
            [obj.WantScattergrams, varargin] = ExtractNamei({'Scattergram','Scattergrams'},varargin);
            [WantAdjust, varargin] = ExtractNamei('Adjust',varargin);
            if WantAdjust
                DefaultNSteps = 200;
                [obj.NSteps, ~] = ExtractNameVali({'NSteps','NStepsApprox'},DefaultNSteps,varargin,'x>10');
                obj.TargetCorrs = obj.RhoControllers;
                FindRhoControllerMatrix(obj);
            end
            % Make sure it is a possible/legal correlation matrix:
            assert(LegalCorrMatrix(obj.RhoControllers),'Error!  The requested correlation matrix is impossible.!');
        end
        
        function Rho = RhoFromController(obj,iRV,jRV,ProposedRhoController)
            % Find the value of Rho that would be obtained between RVs(iRV) and RVs(jRV) when using
            % the proposed RhoController.
            
            % X1 & X2 will hold the generated pairs.
            X1 = repmat(obj.RVPctiles(iRV,:)',obj.NSteps,1);  % X1 is a single column of obj.NSteps^2 X1 values (repeating blocks of obj.NSteps).
            X2 = zeros(obj.NSteps^2,1);
            
            RhoFac = sqrt( 1 - ProposedRhoController^2 );
            
            for ix1=1:obj.NSteps
                MeanAdj = RhoFac * obj.ZPctiles(ix1);
                ThisZ2 = ProposedRhoController * obj.ZPctiles + MeanAdj;
                ThisP2 = normcdf(ThisZ2);
                X2((ix1-1)*obj.NSteps+1:ix1*obj.NSteps) = obj.RVs{jRV}.InverseCDF(ThisP2);
            end
            
            Rho = corr(X1,X2);
            
        end
        
        function FindRhoLimits(obj)
            obj.RhoMin = zeros(obj.NRVs);
            obj.RhoMax = zeros(obj.NRVs);
            for iRV=1:obj.NRVs-1
                X1 = obj.RVPctiles(iRV,:);
                for jRV=iRV+1:obj.NRVs
                    X2 = obj.RVPctiles(jRV,:);
                    obj.RhoMax(iRV,jRV) = corr(X1',X2');
                    X2 = obj.RVPctiles(jRV,end:-1:1);
                    obj.RhoMin(iRV,jRV) = corr(X1',X2');
                end
            end
        end
        
        function RhoController = Find1RhoController(obj,iRV,jRV)
            % Find the value of RhoController for a bivariate normal such that
            % the correlation of percentile-matched RV{iRV} and RV{jRV} (approximated by obj.NSteps)
            % will match the desired TargetCorr
            TargetCorr = obj.TargetCorrs(iRV,jRV);
            ThisMin = obj.RhoMin(iRV,jRV);
            ThisMax = obj.RhoMax(iRV,jRV);
            assert(ThisMin<=TargetCorr && TargetCorr<=ThisMax, ...
                ['Error: cannot attain the requested correlation of ' num2str(TargetCorr) ...
                ' between the RVs ' obj.RVs{iRV}.StringName ' and ' obj.RVs{jRV}.StringName '. ' ...
                'This correlation must be between ' num2str(ThisMin) ' and ' num2str(ThisMax)]);
            if TargetCorr > 0
                LowerBound = 0;
                UpperBound = ThisMax;
            elseif TargetCorr < 0
                LowerBound = ThisMin;
                UpperBound = 0;
            else
                RhoController = 0;
                return;
            end
            fun = @(x) (obj.RhoFromController(iRV,jRV,x)-TargetCorr)^2;
            RhoController = fminbnd(fun,LowerBound,UpperBound,optimset('TolX',.0001,'TolFun',.0001));
        end % Find1RhoController
        
        function FindRhoControllerMatrix(obj)
            % Find the matrix of pairwise RhoController values to obtain the desired TargetCorrs
            % among an arbitrary set of RVs.
            % RVs is a cell array of the random variables.
            % TargetCorrs is an upper-diagonal matrix of the desired correlations.
            % RhoControllerMatrix is an upper-diagonal matrix of the RhoController values
            % required to obtain the desired correlations.
            obj.Pctiles = (1:2:2*obj.NSteps-1) / (2*obj.NSteps);
            obj.ZPctiles = norminv(obj.Pctiles);
            obj.RVPctiles = zeros(obj.NRVs,obj.NSteps);
            for iRV = 1:obj.NRVs
                obj.RVPctiles(iRV,:) = obj.RVs{iRV}.InverseCDF(obj.Pctiles);
            end
            FindRhoLimits(obj);
            obj.RhoControllers = zeros(obj.NRVs);
            for iRV = 1:obj.NRVs-1
                for jRV = iRV+1:obj.NRVs
                    obj.RhoControllers(iRV,jRV) = Find1RhoController(obj,iRV,jRV);
                end
            end
            % copy upper triangular into lower triangular & add ones on diagonal:
            obj.RhoControllers = obj.RhoControllers + triu(obj.RhoControllers,1)' + eye(obj.NRVs);
        end
        
        function rands = Rands(obj,NCases)
            % Generate the correlated random numbers from the RhoController matrix using the normal model.
            mu = zeros(1,obj.NRVs);
            rands = mvnrnd(mu,obj.RhoControllers,NCases);  % NCases rows x NRVs columns; Correlated Normal(0,1)'s.
            rands = normcdf(rands);  % Corresponding CDF values
            for iRV=1:obj.NRVs
                rands(:,iRV) = obj.RVs{iRV}.InverseCDF(rands(:,iRV)); % Arbitrary RVs with the same CDF values.
            end
            if obj.WantHistograms
                obj.Histograms(rands);
            end
            if obj.WantScattergrams
                obj.Scattergrams(rands);
            end
        end
        
        function Histograms(obj,randoms)
            % Display histograms to see whether the marginals look right.
            for iRV=1:obj.NRVs
                figure
                histogram(randoms(:,iRV));
                xlabel(obj.RVs{iRV}.StringName);
            end
        end
        
        function Scattergrams(obj,randoms)
            % Display scattergrams to see that the correlations look right
            obscorrs = corr(randoms);
            for iRV=1:obj.NRVs-1
                for jRV=iRV+1:obj.NRVs
                    figure
                    scatter(randoms(:,iRV),randoms(:,jRV));
                    xlabel(obj.RVs{iRV}.StringName);
                    ylabel(obj.RVs{jRV}.StringName);
                    title(['r = ' num2str(obscorrs(iRV,jRV))]);
                end
            end
        end
        
    end % regular methods
    
end % class
