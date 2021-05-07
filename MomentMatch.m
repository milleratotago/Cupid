function MomentMatch(inDist,outDist,varargin)
    % Adjust the parameters of distribution outDist so that it has the same moment values
    % as the input distribution inDist.  By default, match the first k moments,
    % where k is the number of free parameters in outDist.
    
    % NEWJEFF: Allow for fixed parameters
    % NEWJEFF: Allow for skipped moments by setting their inMoms to NaN
    ParmCodes = outDist.DefaultParmCodes;
    nParms2match = outDist.NDistParms;
    maxMom2compute = nParms2match;
    moms2Compute = 1:maxMom2compute;
    
    inMoms = NaN(maxMom2compute,1);  % these will be observed moments passed to EstMom
    
    for iMom=1:maxMom2compute
        if ismember(iMom,moms2Compute)
            if iMom==1
                inMoms(iMom) = inDist.Mean;
            else
                inMoms(iMom) = inDist.CenMoment(iMom);
            end
        end
    end
    outDist.EstMom(inMoms,ParmCodes);
end

