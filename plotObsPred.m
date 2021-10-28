function [hx, hp] = plotObsPred(X,P,varargin)
    % Plot observed (X) vs predicted (P) distributions.
    % X is a list of data points.
    % P is a Cupid distribution.
    % By default plot PDFs but that can be changed by passing
    %  'CDF' or 'Haz' as varargin
    % hx and hp are the Line structures made by the histogram or plot commands
    
    if (numel(varargin) == 0) ||  strcmpi(varargin{1},'PDF')
        plotType = 1;  % PDF
    elseif strcmpi(varargin{1},'CDF')
        plotType = 2;  % CDF
    elseif strcmpi(varargin{1},'Haz')
        plotType = 3;  % Hazard
    else
        error('Plot type must be PDF, CDF, or Haz');
    end
    
    pCDFs = 0.005:0.01:0.995;
    
    ppts = P.InverseCDF(pCDFs);
    
    switch plotType
        case 1  % PDF
            hx = histogram(X,'normalization','pdf');
            hold on
            ppdf = P.PDF(ppts);
            hp = plot(ppts,ppdf);
        case 2   % CDF
            [xpts, ~, ~, xcdf] = ehaz(X);
            hx = plot(xpts,xcdf);
            hold on
            pcdf = P.CDF(ppts);
            hp = plot(ppts,pcdf);
        case 3   % Hazard
            [xpts, xhaz] = ehaz(X);
            hx = plot(xpts,xhaz);
            hold on
            phaz = P.Hazard(ppts);
            hp = plot(ppts,phaz);
    end  % end switch plotType
    
end
