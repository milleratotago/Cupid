% DemoMLSE

% This demo illustrates how to get standard errors & confidence intervals
% for maximum likelihood parameter estimates using Fisher information.

%% A simple example with the normal distribution:

TrueMu = 0;
TrueSigma = 1;
Nxs = 50;

dist = Normal(TrueMu,TrueSigma);
xs = dist.Random(Nxs,1);

% Fit the dist:
dist.EstML(xs)

% Compute the standard errors of estimation:
[SE, Cov] = dist.MLSE(xs,'rr')

% Compute 95% confidence intervals around the parameter estimates:
zcrit = 1.96;  % for 95% confidence, 2-tailed
mubounds = [dist.mu-zcrit*SE(1), dist.mu+zcrit*SE(1)]
sigmabounds = [dist.sigma-zcrit*SE(2), dist.sigma+zcrit*SE(2)]

%% Repeat this process many times and check how often the bounds capture the estimates:
nSims = 1000;

nMuInbounds = 0;
nSigmaInbounds = 0;

for iSim = 1:nSims
    dist.ResetParms([TrueMu, TrueSigma]);
    xs = dist.Random(Nxs,1);
    dist.EstML(xs);
    [SE, Cov] = dist.MLSE(xs,'rr');
    mubounds = [dist.mu-zcrit*SE(1), dist.mu+zcrit*SE(1)];
    sigmabounds = [dist.sigma-zcrit*SE(2), dist.sigma+zcrit*SE(2)];
    if (TrueMu>mubounds(1)) && (TrueMu<mubounds(2))
        nMuInbounds = nMuInbounds + 1;
    end
    if (TrueSigma>sigmabounds(1)) && (TrueSigma<sigmabounds(2))
        nSigmaInbounds = nSigmaInbounds + 1;
    end
end
PrMuInbounds = nMuInbounds / nSims
PrSigmaInbounds = nSigmaInbounds / nSims