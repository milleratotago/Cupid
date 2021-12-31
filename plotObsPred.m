function [hx, hp] = plotObsPred(X,P,plotType,varargin)
    % Plot observed (X) vs predicted (P) distributions.
    % X is a list of data points.
    % P is a Cupid distribution.
    % By default plot PDFs but that can be changed by passing
    %  'CDF' or 'Haz' as varargin
    % hx and hp are the Line structures made by the histogram or plot commands
    
    astruc = struct();
    if nargin<=2
        plotType = 1;
    else
        if numel(varargin)>0
            astruc = varargin{1};
        end
    end
%     if numel(varargin) == 0
%         plotType = 1;  % PDF
%     elseif isnumeric(varargin{1})
%         plotType = varargin{1};
%     elseif strcmpi(varargin{1},'PDF')
%         plotType = 1;  % PDF
%     elseif strcmpi(varargin{1},'CDF')
%         plotType = 2;  % CDF
%     elseif strcmpi(varargin{1},'Haz')
%         plotType = 3;  % Hazard
%     else
%         error('Plot type must be PDF, CDF, or Haz');
%     end
    
    pCDFs = [0.001 0.005:0.01:0.995 0.999];
    
    ppts = P.InverseCDF(pCDFs);
    
    switch plotType
        case 1  % PDF
            hx = histogram(X,'normalization','pdf');
            hold on
            % Plotting PDFs may not look right if the PDF
            % changes sharply within a bin.
            % ppdf = P.PDF(ppts);
            % hp = plot(ppts,ppdf);
            % For that reason, it is better to plot
            % the probabilities in the different bins.
            BinCDFs = P.CDF(hx.BinEdges);  % CDFs at the top of each bin
            BinPDFs = diff(BinCDFs) / hx.BinWidth;
            NBins = length(hx.BinEdges) - 1;
            % Now represent each bin with 2 X values
            % for its top and bottom edges, both having same PDF
            % to get a flat line.
            fnlBinEdges = zeros(2*NBins-1,1);
            fnlBinPDFs = zeros(2*NBins-1,1);
            fnliBin = 0;
            for iBin=1:NBins
                fnliBin = fnliBin + 1;
                fnlBinEdges(fnliBin) = hx.BinEdges(iBin);
                fnlBinPDFs(fnliBin) = BinPDFs(iBin);
                fnliBin = fnliBin + 1;
                fnlBinEdges(fnliBin) = hx.BinEdges(iBin+1) - eps;
                fnlBinPDFs(fnliBin) = BinPDFs(iBin);
           end
            plot(fnlBinEdges,fnlBinPDFs);
        case 2   % CDF
            hx = histogram(X,'normalization','cdf');
%             [xpts, ~, ~, xcdf] = ehaz(X);
%             hx = plot(xpts,xcdf);
            hold on
            pcdf = P.CDF(ppts);
            hp = plot(ppts,pcdf);
        case 3   % Hazard
            [xpts, xhaz] = ehaz(X);
            hx = plot(xpts,xhaz);
            hold on
            phaz = P.Hazard(ppts);
            hp = plot(ppts,phaz);
        otherwise
            error('Unrecognized plotType');
    end  % end switch plotType
%     for iFun=1:numel(afun)
%         afun{iFun}();
%     end
    sFields = fieldnames(astruc);
    nFields = numel(sFields);
    for iField = 1:nFields
        func = str2func(sFields{iField});
        func(astruc.(sFields{iField}));
    end
%     if isfield(astruc,'legend')
%         legend('observed','predicted');
%     end
end
