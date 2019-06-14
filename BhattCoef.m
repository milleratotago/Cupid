function c = BhattCoef(dist1,dist2)
    % Bhattacharyya coefficient measuring similarity of two probability distributions.
    % See e.g. https://en.wikipedia.org/wiki/Bhattacharyya_distance
    if (dist1.DistType == 'd') && (dist2.DistType == 'd')
        commonXs = intersect(dist1.DiscreteX,dist2.DiscreteX);
        prodPDFs = dist1.PDF(commonXs) .* dist2.PDF(commonXs);
        c = sum(sqrt(prodPDFs));
    elseif (dist1.DistType == 'c') && (dist2.DistType == 'c')
        fn2int = @(x)(sqrt(dist1.PDF(x).*dist2.PDF(x)));
        minx = max(dist1.LowerBound,dist2.LowerBound);
        maxx = min(dist1.UpperBound,dist2.UpperBound);
        c = integral(fn2int,minx,maxx);
    else
        error('Bhattacharyya coefficient can only be computed for 2 discrete or 2 continuous distributions');
    end
end
