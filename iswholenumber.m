function thisval=iswholenumber(X)
thisval = isreal(X) && (rem(X,1)<1e-12);
% if ~thisval
%     fprintf('%f is not a whole number\n',abs(round(X)-X));
% end
end
