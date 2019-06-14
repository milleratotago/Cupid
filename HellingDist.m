function d = HellingDist(dist1,dist2)
    % Hellinger distance between two probability distributions.
    % See e.g. https://en.wikipedia.org/wiki/Bhattacharyya_distance
    d = sqrt(1 - BhattCoef(dist1,dist2) );
end
