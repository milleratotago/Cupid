% Examples of Probit analysis


% The first part of this demo illustrates analyses for the yes/no task,
% and the second part illustrates mAFC.


%% ********************** Yes/no task

% Generate some artificial data to be used in illustrating the YN analyses.

mydist = Normal(50,5);                 % Hypothesized underlying distribution
Cs = mydist.InverseCDF(.10:.20:.90);   % Values at which to run tests
Ns = 100*ones(size(Cs));               % Repeat the test 100 times at each value
pSmaller = [.05 .31 .51 .71 .95];      % Approximately the CDF values, with some error thrown in so that the
                                       % distribution parameters will need some adjustment
pGreater = 1 - pSmaller;               % Proportion of times an X value is greater than Cs(i)
Gs = Ns .* pGreater;                   % Number of times observed value is greater than Cs(i)

%% Illustrate the command YNProbitLnLikelihood:
%  Compute the log-likelihood of a dataset (in this case, an artificial dataset)
%  under the hypothesized distribution.
mydist.YNProbitLnLikelihood(Cs,Ns,Gs)


%% Illustrate the command EstProbitYNML

%  Estimate the mean and sd of the best-fitting normal
%  based on these artifical data, by maximizing likelihood.
mydist.EstProbitYNML(Cs,Ns,Gs)

mydist.YNProbitLnLikelihood(Cs,Ns,Gs)   % Note that the likelihood is a little larger with the new parameters


mydist = Normal(50,5);  % Reset to the starting distribution.

%% Illustrate the command YNProbitChiSq:
%  Compute the chi-square goodness of fit of the artificial data
%  under the distribution.
mydist.YNProbitChiSq(Cs,Ns,Gs)

%% Illustrate the command EstProbitYNChiSq
%  Estimate the mean and sd of the best-fitting normal
%  based on these artifical data, by minimizing chi-square.
mydist.EstProbitYNChiSq(Cs,Ns,Gs)

mydist.YNProbitChiSq(Cs,Ns,Gs)   % Note that the chi-square is a little smaller with the new parameters



%% ********************** mAFC task

% Generate some artificial data to be used in illustrating the mAFC analyses.

mAFC = 3;  % 3-alternative forced-choice task

mydist = Normal(50,5);                   % Hypothesized underlying distribution
pCorrect = pSmaller + (1-pSmaller)/mAFC; % pCorrect = pSmaller plus correct guesses
Gs = Ns .* pCorrect;                     % Number of times the observed response is correct

%% Illustrate the command mAFCProbitLnLikelihood:
%  Compute the log-likelihood of the artificial data
%  under the distribution.
mydist.mAFCProbitLnLikelihood(mAFC,Cs,Ns,Gs)

%% Illustrate the command EstProbitmAFCML
%  Estimate the mean and sd of the best-fitting normal
%  based on these artifical data, by maximizing likelihood.
mydist.EstProbitmAFCML(mAFC,Cs,Ns,Gs)

mydist.mAFCProbitLnLikelihood(mAFC,Cs,Ns,Gs)   % Note that the likelihood is a little larger with the new parameters

mydist = Normal(50,5);  % Reset to the starting distribution.

%% Illustrate the command mAFCProbitChiSq:
%  Compute the log-likelihood of the artificial data
%  under the distribution.
mydist.mAFCProbitChiSq(mAFC,Cs,Ns,Gs)


%% Illustrate the command EstProbitmAFCChiSq
%  Estimate the mean and sd of the best-fitting normal
%  based on these artifical data, by minimizing chi-square.
mydist.EstProbitmAFCChiSq(mAFC,Cs,Ns,Gs)

mydist.mAFCProbitChiSq(mAFC,Cs,Ns,Gs)   % Note that the chi-square is a little smaller with the new parameters


