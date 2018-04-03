function thispr = PrXGTY(xDist,yDist)
% Compute prob(X>Y) for any two independent random variables. LIMITED AT THIS POINT TO CONTINUOUS RVS.
thispr = integral(@(x) PDF(xDist,x).*CDF(yDist,x),xDist.LowerBound,xDist.UpperBound);
end
