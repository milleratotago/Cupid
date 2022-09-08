%% Demo of Probit Analysis
% This demo shows how to use Cupid for "probit analysis", which
% I describe using the terminology of psychophysical experiments.
% For further information, look in the PDF documentation file.
%
% The first part of the demo illustrates analyses for the yes/no task,
% and the second part illustrates m-alternative forced-choice (mAFC) tasks.


%% Yes/no task, maximum-likelihood fit
%
% Generate some artificial data for illustrating the YN analyses.
% Normally you would have your own data, which should be in an analogous format.

mydist = Normal(50,5);                 % Hypothesized underlying distribution of X's
Cs = mydist.InverseCDF(.10:.20:.90);   % A set of stimulus values at which to run tests.
                                       % At each test the result will be "Yes" if Cs>X and "no" if Cs<X.
Ns = 100*ones(size(Cs));               % Repeat the test 100 times at each value
PrCgreaterthanX = ...
    [.05 .31 .51 .71 .95];             % The proportion of times that X is less than Cs(i) for each i.
                                       % These are approximately the CDF values of the Cs, but not exactly so that the
                                       % distribution parameters will need some adjustment
Gs = Ns .* PrCgreaterthanX;            % Number of observations of Cs(i) greater than X

%%
% The sample data have now been generated.
%%
% Next we try to fit the sample data with the assumed underlying Normal distribution.
% The fitting process estimates the mean and sd of the normal based on these artifical data, by maximizing likelihood.
mydist.EstProbitYNML(Cs,Ns,Gs)
%%
% If you want, you can directly compute the likelihood of the sample data under the fitted distribution.
mydist.YNProbitLnLikelihood(Cs,Ns,Gs)

%% Yes/no task, minimum-chi-square fit
% Alternatively, it is possible to estimate parameters by minimizing chi-square.
mydist.EstProbitYNChiSq(Cs,Ns,Gs)


%% Spearman-Karber analysis.
% We will re-use the Cs and PrCgreaterthanX values defined above, but this
% time we will generate observed counts using the binomial distribution.
NTrialsPerC = 100;
CgreaterCount = binornd(NTrialsPerC,PrCgreaterthanX);
ClesserCount = NTrialsPerC - CgreaterCount;
monoCDFs = SpearKarDist.monotonize(ClesserCount,CgreaterCount);

% Augment Cs and monoCDFs, if necessary, so that monoCDFs includes 0 and 1.
% The choice of corresponding min and max Cs is somewhat arbitrary.
[Cs,monoCDFs] = SpearKarDist.Stretch01(Cs,monoCDFs);

% Create a SpearKar distribution based on these values
SK = SpearKarDist(Cs,monoCDFs);

% Compute nonparametric estimates of the mean, variance, etc
SK.Mean
SK.SD
SK.Variance
SK.Median

%% mAFC task, maximum-likelihood fit
% Generate some artificial data to be used in illustrating the mAFC analyses.
% In these tasks each response is "correct" or "incorrect".
% Each response is correct if (Cs>X) || ( (Cs<X) && (there is a correct guess from mAFC alternatives) )

mAFC = 3;  % 3-alternative forced-choice task

mydist = Normal(50,5);            % Hypothesized underlying distribution
pCorrect = PrCgreaterthanX + ...
    (1-PrCgreaterthanX)/mAFC;     % pCorrect = PrCgreaterthanX plus correct guesses
Gs = Ns .* pCorrect;              % Number of times the observed response is correct

%%
% Now we can use these simulated data to estimate the mean and sd of the best-fitting normal
% based on these artifical data, by maximizing likelihood.
mydist.EstProbitmAFCML(mAFC,Cs,Ns,Gs)


%% mAFC task, minimum-chi-square fit
% Alternatively, it is possible to estimate parameters by minimizing chi-square.
mydist.EstProbitmAFCChiSq(mAFC,Cs,Ns,Gs)


