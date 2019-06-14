function d = BhattDist(dist1,dist2)
    % Bhattacharyya distance between two probability distributions.
    % See e.g. https://en.wikipedia.org/wiki/Bhattacharyya_distance
    d = -log(BhattCoef(dist1,dist2));
end
