classdef ExConvFFTc < ConvolveFFTc
    % ExConvFFTc(BasisRV1,exp_mean)
    % Convolution via ConvolveFFTc of any continuous BasisRV1 with an ExponenMn

    methods
        
        function obj=ExConvFFTc(BasisRV1, exp_mean,varargin)
            % StartParmsMLE assumes BasisRV1 is first.
            obj=obj@ConvolveFFTc(BasisRV1,ExponenMn(exp_mean),varargin{:});
            obj.ParmTypes = [BasisRV1.ParmTypes 'r'];
            obj.DefaultParmCodes = [BasisRV1.DefaultParmCodes 'r'];
            obj.NDistParms = BasisRV1.NDistParms + 1;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            for iParm = 1:BasisRV1.NDistParms
                obj.ParmNames{iParm} = BasisRV1.ParmNames{iParm};
            end
            obj.ParmNames{obj.NDistParms} = 'exmean';
        end
        
       function parms = StartParmsMLE(obj,X)

            % Start with exponential mean suggested by upper tail of X's.
            nX = numel(X);
            sortX = sort(X);
            hiP = (nX - 0.5) / nX;
            hiT = sortX(end);
            tailLen = ceil(0.05*nX);
            if tailLen > 10
                tailLen = 10;   % Don't use a tail longer than 10 X's.
            elseif tailLen < 2
                tailLen = 2;    % tail must include at least 2 X's.
            end
            lesserT = sortX(end-(tailLen-1));
            lesserP = (nX - tailLen + 0.5) / nX;
            % Next line derived from difference in two CDF values for an exponential
            estRate = (log(1-lesserP) - log(1-hiP)) / (hiT - lesserT);
            estExpMean = 1/estRate;  %  * haztp1 / haztp2;
            % Note: I think it should be possible to improve the estimate of Rate by looking
            % at the CDF of BasisRV1 for hiT and lesserT but I have not yet worked out how.
            
            % Now we can estimate the remaining mean & sd of the non-exponential
            % A problem can arise here if there are large outliers--namely,
            % nonExpMean and/or nonExpVar can sometimes be estimated to be negative.
            % Here we make sure nonExpVar is positive; individual distributions
            % will have to take care of the mean if a negative is not allowed.
            nonExpMean = mean(X) - estExpMean;
            nonExpVar = max([eps, 0.1*nonExpMean, std(X)^2 - estExpMean^2]);
            % CLUGE: Temporarily increase tolerance of parameter searches since we are just looking for starting values:
            IncreaseFactor = 1e4;
            obj.BasisRV{1}.SearchOptions.TolFun = obj.BasisRV{1}.SearchOptions.TolFun * IncreaseFactor;
            obj.BasisRV{1}.SearchOptions.TolX = obj.BasisRV{1}.SearchOptions.TolX * IncreaseFactor;
            obj.BasisRV{1}.EstMom([nonExpMean, nonExpVar]);
            obj.BasisRV{1}.SearchOptions.TolFun = obj.BasisRV{1}.SearchOptions.TolFun / IncreaseFactor;
            obj.BasisRV{1}.SearchOptions.TolX = obj.BasisRV{1}.SearchOptions.TolX / IncreaseFactor;
            parms = [obj.BasisRV{1}.ParmValues estExpMean];  % 10 is mean of fast exponential
       end
       
    end % methods

end % classdef
