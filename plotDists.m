function [x, y] = plotDists(Dists,PdfCdf,varargin)   % NEWJEFF: Undocumented
    % Utility function to plot one PDF or CDF for each of several Cupid probability distributions.
    % Dists: a cell array of the probability distributions to be plotted.
    % PdfCdf = 1/2 to plot PDF or CDF
    % Optional arguments:
    %  nPoints: number of x values to compute/plot for each distribution
    %  The remainder of varargin is a list of parameters to be passed
    %  to the plot command for control of line type, thickness, etc.
    %  The first 1/nDists of these parameters are passed when plotting the first distribution,
    %  the next 1/nDists are passed when plotting the 2nd distribution, ...
    [nPoints, varargin] = ExtractNameVali('nPoints',100,varargin);
    nDists = numel(Dists);
    x = nan(nDists,nPoints);
    y = nan(nDists,nPoints);
    ttlNPlotParms = numel(varargin);
    nPlotParmsPerDist = floor(ttlNPlotParms/nDists);
    distNames = {};
    hold on;
    for iDist=1:nDists
        if nPlotParmsPerDist>0
            plotParms = varargin((iDist-1)*nPlotParmsPerDist+1:iDist*nPlotParmsPerDist);
        else
            plotParms = {};
        end
        thisDist = Dists{iDist};
        x(iDist,:) = linspace(thisDist.LowerBound,thisDist.UpperBound,nPoints);
        if PdfCdf == 1
            y(iDist,:) = thisDist.PDF(x(iDist,:));
        elseif PdfCdf == 2
            y(iDist,:) = thisDist.CDF(x(iDist,:));
        end
        plot(x(iDist,:),y(iDist,:),plotParms{:});
        distNames = [distNames {thisDist.StringName}];
    end
    legend(distNames);
end
