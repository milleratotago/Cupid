% Probit stuff to be added

Cs = .05:.1:.95;    % These are RV values within the range of mydist
Ns = 100*ones(size(Cs));  % Integers
NGs = 105 - (1:10).*Ns/10;
mydist.YNProbitLnLikelihood(Cs,Ns,NGs)
NGs2 = 100 - NGs;
mydist.YNProbitLnLikelihood(Cs,Ns,NGs2)

mydist.mAFCProbitLnLikelihood(m,N,Cs,Ns,NGs)
