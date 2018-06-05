function thispr = PrXGTY(xDist,yDist)
    % Compute prob(X>Y) for any two independent random variables.
    % In this version, at least one RV must be continuous so that
    % there is no possibility of ties.
    
    if (xDist.DistType == 'c') && (yDist.DistType == 'c')
        thispr = integral(@(x) PDF(xDist,x).*CDF(yDist,x),xDist.LowerBound,xDist.UpperBound);
    elseif (xDist.DistType == 'd') && (yDist.DistType == 'c')
        thispr = sumofintegrals(xDist,yDist);
    elseif (xDist.DistType == 'c') && (yDist.DistType == 'd')
        thispr = 1 - sumofintegrals(yDist,xDist);
    else
        error('PrXGTY cannot handle two discrete distributions.');
        % NWJEFF: Could at least handle two discrete distributions if there were no ties.
    end
end

function thissum = sumofintegrals(xDiscrete,yContinuous)
    % Compute the probability that a discrete random variable X
    % is larger than a continuous random variable Y.
    yCDF = yContinuous.CDF(xDiscrete.DiscreteX);
    thissum = sum( xDiscrete.DiscretePDF .* yCDF);
%     return
%     thissum = 0;
%     for iX=1:xDiscrete.NValues
%         thissum = thissum + xDiscrete.DiscretePDF(iX) * yContinuous.CDF(xDiscrete.DiscreteX(iX));
%     end
end
