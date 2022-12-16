classdef FrechetPos < Frechet
    % Version of the Frechet where the minimum value must be >= 0

    properties
        Property1
    end

    methods

        function obj=FrechetPos(varargin)
            obj=obj@Frechet(varargin{:});
            obj.FamilyName = 'FrechetPos';
            obj.NDistParms = 3;
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.CDFNearlyOne = 0.999999;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            obj.fsolveoptions = optimset('Display','off');
            obj.maxCDF = 0.9999;
            obj.maxUpperBound = 50000;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{1},varargin{2},varargin{3}]);
                otherwise
                    ME = MException('FrechetPos:Constructor', ...
                        'FrechetPos constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end

        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.shape = newparmvalues(1);
            obj.scale = newparmvalues(2);
            obj.minval = abs(newparmvalues(3));
            ReInit(obj);
        end
        
        function parms = StartParmsMLE(~,X)
            % Guess the min value as a litle less than min(X), >=0
            Xmin = min(X);
            Xsd = std(X);
            MinGuess = max(0,Xmin - Xsd/(100*numel(X)));
            % Subtract the min so that the distribution starts at 0.
            Y = X - MinGuess;
            obsmedian = median(Y);
            estscale = obsmedian;
            obs90pct = prctile(Y,90);
            ErrFn = @(x) 0.90 - Frechet.frechcdf(obs90pct,x,estscale,0);
            estshape = fzero(ErrFn,3);
            parms = [estshape, estscale, MinGuess];
        end

    end % methods

end % classdef