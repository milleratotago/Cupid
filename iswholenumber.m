function thisval = iswholenumber(X)
    % Sometimes rem is almost 0 or 1
    thisval = isreal(X) && ( (rem(X,1)<1e-10) || (rem(X,1)>1-1e-10) );
end
