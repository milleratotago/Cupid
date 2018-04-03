% Script to give a brief demonstration of what CUPID does & how to use it (approx 100 lines).
% It is best to step through it with successive clicks on "Run Section",
% checking the commands, comments, and output of each section one by one.

%% Create a distribution object by indicating its name and parameters, like this:
IQdist = Normal(100,10);
WaitingTime = Exponential(.3);
Cost = Weibull(3,10,2);
% The documentation lists dozens of basic distributions that have been implemented so far.
% It is easy to add new distributions for which you know at least either the PDF or the CDF.


%% These new variables are statistical distributions, and you can retrieve their properties, e.g.
IQmean = IQdist.Mean
IQpctileof120 = IQdist.CDF(120)
IQat99thpctile = IQdist.InverseCDF(.99)
WaitingTimeVar = WaitingTime.Variance;
CostMedian = Cost.Median;
% The documentation lists the many functions that can be computed for each distribution.

% Because of the object-oriented approach, all of these properties can be computed
% for any new distribution, even if only the PDF or CDF is available.

% Here are some more examples of supported functions.

%% Request/plot PDF, CDF.
IQpdf = IQdist.PDF(50:150);
plot(50:150,IQpdf);
IQdist.PlotDens;  % PlotDens is a command to make plots of the PDF and CDF of any distribution.

%% Parameter estimation.
% The parameters of any distribution can be estimated according to various criteria,
% including: maximum likelihood, method of moments, chi-square goodness of fit,
% matching certain desired percentiles, etc.
% As one example, here is the command to estimate the parameters of the IQdist
% distribution so that the IQdist scores of 80 and 120 are at the 40th and 60th
% percentiles, respectively.
IQdist.EstPctile([80 120],[0.4 0.6])
fprintf('This gives a new distribution with mean and sd of %f and %f,\n',IQdist.Mean,IQdist.SD)
fprintf(' for which scores 80 and 120 have CDF of %f and %f as requested.\n\n',IQdist.CDF(80),IQdist.CDF(120));

%% During parameter estimation, it is also possible to hold one or more parameter values fixed
% by including an optional "ParmCodes" parameter.
% For example, suppose we want to hold the mean fixed at 100 and adjust the SD
% so that an IQdist of 140 is at the 99th percentile:
IQdist = Normal(100,10);
IQdist.EstPctile([140],[0.99],'fr')  % 'fr' says to fix the first parameter and vary the second parameter as a real number.
fprintf('This gives a new distribution with mean and sd of %f and %f,\n',IQdist.Mean,IQdist.SD)
fprintf(' for which 140 has a CDF of %f as requested.\n\n',IQdist.CDF(140));


%% Derived (transformed) distributions, 1.
% In addition to having many standard distributions, CUPID has many distribution classes
% that represent new distributions formed as a transformation of some other distributions.
% Examples:
% Here is the distribution of the square-root of a gamma distribution:
SqrtG = SqrtTrans(RNGamma(3,1));
% A transformed distribution can be treated just like any other, e.g.:
[ SqrtG.Mean SqrtG.Variance]
SqrtG.EstMom([2.3 .09])
[ SqrtG.Mean SqrtG.Variance]

% As another example, here is the distribution of the log of the previous Weibull distribution:
LogCost = LogTrans(Cost);
[LogCost.Mean LogCost.SD]

% The documentation lists many supported transformations.

%% Other derived distributions, 2.
% Other CUPID classes support distributions that are derived from multiple
% underlying distributions, including distributions of sums (convolutions),
% differences, mixtures, order statistics, truncated distributions, ...
% For example, here is the distribution of the minimum (1st order statistic) in a sample
% of 3 independent scores, one of which is normal, one of which is exponential,
% and one of which is uniform.
Min1 = Order(1,Normal(0,1),Exponential(.3),Uniform(0,1));
Min1.Median
Min1.PlotDens

%% There is of course a short-cut for order statistics from an IID sample:
% For example, "Ord2of8" is defined as the second order statistic
% from a sample of 8 Uniform(0,1) random variables:
Ord2of8 = OrderIID(2,8,Uniform(0,1))
fprintf('This gives a new distribution with mean and sd of %f and %f,\n',Ord2of8.Mean,Ord2of8.SD)
Ord2of8.PlotDens

%% Because of the OO approach, any derived distribution is itself a new distribution,
% and so a further new distribution may be derived from it:
TruncOrdP = TruncatedX(Ord2of8,0,.05)
fprintf('This gives a new distribution with mean and sd of %f and %f,\n',TruncOrdP.Mean,TruncOrdP.SD)
TruncOrdP.PlotDens

%% As the next example shows, these distribution specifications can be composed, as in f(g(x)).
% In this example we start with the Uniform(0,1) distribution, truncate it so that it is
% restricted to the range p<.05, and then take the 2nd order statistic of 8 such RVs.
% Do you think that gives the same distribution as TruncOrdP formed earlier?
OrdTruncP = OrderIID(2,8,TruncatedX(Uniform(0,1),0,.05))
fprintf('This gives a new distribution with mean and sd of %f and %f,\n',OrdTruncP.Mean,OrdTruncP.SD)
OrdTruncP.PlotDens
% Answer: no, this is quite a different distribution!

