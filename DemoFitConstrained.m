%% Fitting distributions with constraints
% This demo script gives some examples of how to fit one or more distributions
% with certain contraints on the distribution(s) parameter(s).
%
% This type of fitting requires you to write special functions to enforce the constraints,
% as will be illustrated with numerous examples.  The example constraint
% functions are all in the file DemoFitConstrainedFns.m.
%
% **NOTE**: In each example, the first step is to generate a fake dataset
% so that there are data available for fitting. Presumably you would already
% have the data that you wanted to fit, so you would skip this step.

%% Example: Maximum likelihood fit of a normal distribution with the constraint sigma>0.1*mu.
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=100;  % arbitrary sample size
Datasets{1} = Normal(55,4.5).Random(n1,1);
%%
% Note that the constraint is not truly satisfied because the true sigma is less than 5.5=0.1*mu.
%
% Set some parameters defining the desired fit:
Dists{1} = Normal(45,4);       % The distribution to be fit.  These parameter values are not actually used.
sErrorFn = '-LnLikelihood';    % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [40, sqrt(5)];  % Starting values for the two parameters expected by DemoFitConstraintFn0.
ConstraintFn = @DemoFitConstrainedFns.Fn0;
%%
% The user-supplied constraint function converts fminsearch's suggested parameter values
% into suggested parameter values for the distributions that are being fit, implementing the constraint.
% Look at the example functions to see exactly how that is done.

%%
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distribution is ' Dists{1}.StringName]);

%%
% What would the fit have been without the constraint?
nocon = Dists{1}.EstML(Datasets{1});
disp(['In comparison, without the constraint the estimated distribution is:' nocon]);


%% Example: Maximum likelihood fit of two normal distributions with the constraint that they have the same sigma.
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=100; n2=90;  % arbitrary sample sizes
Datasets{1} = Normal(10,10).Random(n1,1);
Datasets{2} = Normal(20,12).Random(n2,1);
%%
% Note that the constraint is not truly satisfied because the true sigmas are not equal.
%
% Set some parameters defining the desired fit:
Dists{1} = Normal(15,11);  % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = Normal(15,11);
sErrorFn = '-LnLikelihood';  % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [15, 15, sqrt(11)];      % Starting values for the three parameters expected by DemoFitConstraintFn1.
ConstraintFn = @DemoFitConstrainedFns.Fn1;
%%
% We supply a different constraint function to implement these constraints.
%
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);


%% Example: Maximum likelihood fit of three distributions with equality across distributions of certain parameters.
% For this example, we fit three distributions: a normal, an exponential, and an ex-Gaussian, with the
% constraint that the parameters of the ex-Gaussian are equal to those of the normal and exponential.
%
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=100; n2=90; n3 = 130;  % arbitrary sample sizes
Datasets{1} = Normal(300,30).Random(n1,1);
Datasets{2} = Exponential(1/100).Random(n2,1);
Datasets{3} = Normal(320,32).Random(n3,1) + Exponential(1/100).Random(n3,1);
%%
% Note that the constraint is not truly satisfied because the parameters of the exGaussian
% are not identical to those of the normal and exponential.
%
% Set some parameters defining the desired fit:
Dists{1} = Normal(15,11);       % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = ExponenMn(100);      % Use the exponential whose parameter is the mean, not rate.
Dists{3} = ExGauMn(15,11,100);  % Use the exGaussian whose exponential parameter is the mean, not rate.
sErrorFn = '-LnLikelihood';     % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [290, sqrt(23), sqrt(100)];  % Starting values for the three parameters expected by DemoFitConstraintFn2.
ConstraintFn = @DemoFitConstrainedFns.Fn2;
%%
% Yet another constraint function.
%
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ', ' Dists{2}.StringName ', and ' Dists{3}.StringName]);



%% Example: ChiSquare fit of two normal distributions with the constraint that they have the same sigma
% This example requires construction of more complicated Datasets, because the GofFChiSq error function
% requires both bin specifications and counts as input data.
%
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=400; n2=590;  % arbitrary sample sizes
Scores1 = Normal(10,10).Random(n1,1);
Scores2 = Normal(20,12).Random(n2,1);
%%
% Note that the constraint is not truly satisfied because the true sigmas are not equal.
%
% Use histcounts to bin and count the scores, since GofFChiSq works
% with bins and probabilities, not raw scores.
[counts1,bins1] = histcounts(Scores1);
[counts2,bins2] = histcounts(Scores2);
%%
% Make the datasets that will be passed to GofFChiSq.
% Each dataset specifies GofFChiSq's parameters of BinUpperBounds and BinProbs, in that order.
Datasets{1}{1} = bins1(2:end);
Datasets{1}{2} = counts1/sum(counts1);
Datasets{2}{1} = bins2(2:end);
Datasets{2}{2} = counts2/sum(counts2);
%%
% Now we are done generating the fake datasets.
%
% Set some parameters defining the desired fit:
Dists{1} = Normal(15,11);  % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = Normal(15,11);
sErrorFn = 'GofFChiSq';  % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [15, 15, sqrt(11)];      % Starting values for the three parameters expected by DemoFitConstraintFn1.
ConstraintFn = @DemoFitConstrainedFns.Fn1;

%%
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);



%% Example: Fit two logistic distributions with a constraint on the sum of their CDFs at a certain point.
% In this example we fit two logistic distributions to two 2AFC psychometric functions with
% the constraint described by Bausenhart et al, Behav Res (2012) 44:1157–1174
% DOI 10.3758/s13428-012-0207-z
% The constraint is that the CDFs of these 2 distributions evaluated at s=50 should sum to 1.0
% This example also requires construction of complicated Datasets, because
% the mAFCProbitLnLikelihood error function requires multiple inputs.
%
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
m = 2;       % Use a 2AFC task for this example.
pGuess = 1/m;
s = 50;      % This example uses a constant of s=50 (also used in DemoFitConstraintFn4)
Constants = [10:10:90];     % Stimulus values used in psychometric testing
NTrials = 50*ones(size(Constants));   % Number of tests at each constant
CDF1 = Logistic(45,10).CDF(Constants);
CDF2 = Logistic(55,10).CDF(Constants);
PCorr1 = CDF1 + (1-CDF1)*pGuess;
PCorr2 = CDF2 + (1-CDF2)*pGuess;
NCorrect1 = binornd(NTrials,PCorr1);
NCorrect2 = binornd(NTrials,PCorr2);
%%
% Note: I think that the constraint is truly satisfied in this case.
%
% Now finish constructing the datasets that will be passed to mAFCProbitLnLikelihood.
% Each dataset specifies mAFCProbitLnLikelihood's required parameters, in order.
% Here is the function definition as a reminder of what datasets are required:
% function thisval=mAFCProbitLnLikelihood(obj, PassM, Constants, NTrials, NCorrect)
Datasets{1}{1} = m;  % the value of PassM
Datasets{1}{2} = Constants;
Datasets{1}{3} = NTrials;
Datasets{1}{4} = NCorrect1;
Datasets{2}{1} = m;
Datasets{2}{2} = Constants;
Datasets{2}{3} = NTrials;
Datasets{2}{4} = NCorrect2;
%%
% Now we are done generating the fake datasets.
%
% Set some parameters defining the desired fit:
Dists{1} = Logistic(50,10);  % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = Logistic(50,10);
sErrorFn = '-mAFCProbitLnLikelihood';   % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [15, 15, sqrt(11)];      % Starting values for the three parameters expected by DemoFitConstraintFn4.
ConstraintFn = @DemoFitConstrainedFns.Fn4;

%%
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);
disp('The constraint implies that the CDFs of these distributions evaluated at s=50 should sum to 1.0:');
disp(['The actual sum is ' num2str(Dists{1}.CDF(s) + Dists{2}.CDF(s))]);


%% Example: Constrain properties of the distributions rather than parameters.
% Example: Constrain some properties of the distributions--here, their means and SDs--
% rather than constraining their parameters.
% In this example, fminsearch can suggest any real parameter values that it wants, and this
% could cause the search to fail if an illegal parameter value is suggested (e.g., negative variance)
% See the following example for a more bullet-proof approach.
%
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=200; n2=190;  % arbitrary sample sizes
TrueDist1 = Laplace(105,21);
fprintf('The true mean and standard deviation of Dist1 are %f and %f.\n',TrueDist1.Mean,TrueDist1.SD);
TrueDist2 = RNGamma(10,.10);   % Second true distribution is Gamma
fprintf('The true mean and standard deviation of Dist2 are %f and %f.\n',TrueDist2.Mean,TrueDist2.SD);
Datasets{1} = TrueDist1.Random(n1,1);   % Data from first true distribution.
Datasets{2} = TrueDist2.Random(n2,1);   % Data from first true distribution.
%%
% Note that the constraint is not truly satisfied because the distributions do not
% have identical means & std devs.
% Now we are done generating the fake datasets.
%
% Set some parameters defining the desired fit:
Dists{1} = Laplace; % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = RNGamma;
sErrorFn = '-LnLikelihood';   % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [100,20,10,.10];      % Starting values for the four parameters expected by DemoFitConstraintFn5.
ConstraintFn = @DemoFitConstrainedFns.Fn5;

%%
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);
disp('The constraint implies that the fitted distributions should have equal means and std devs:');
fprintf('The estimated mean and standard deviation of Dist1 are %f and %f.\n',Dists{1}.Mean,Dists{1}.SD);
fprintf('The estimated mean and standard deviation of Dist2 are %f and %f.\n',Dists{2}.Mean,Dists{2}.SD);

%% Example: Constrain properties of the distributions rather than parameters (more bullet-proof).
% Example: Constrain two distributions to produce the same mean and standard deviation.
% This example relies on the distribution-specific functions ParmsToReals and RealsToParms
% to constrain each distribution''s parameter values to their possible ranges.
% This requires a bit more elaborate programming, but it is better than the previous example
% because it ensures that no illegal parameter combinations are tried.
%
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=200; n2=190;  % arbitrary sample sizes
TrueDist1 = Normal(105,18);
TrueDist2 = RNGamma(11,1/8);
Datasets{1} = TrueDist1.Random(n1,1);   % Data from distributions.
Datasets{2} = TrueDist2.Random(n2,1);
%%
% Note that the constraint is not truly satisfied because these distributions
% do not have the same variance.
% Now we are done generating the fake datasets.
%
% Set some parameters defining the desired fit:
Dists{1} = Normal(100,10);
Dists{2} = RNGamma(10,1/10);
sErrorFn = '-LnLikelihood';   % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [Dists{1}.ParmsToReals(Dists{1}.ParmValues) Dists{2}.ParmsToReals(Dists{2}.ParmValues)];      % Starting values for the parameters expected by DemoFitConstraintFn7.
ConstraintFn = @DemoFitConstrainedFns.Fn7;

%%
% Now perform the actual fit:
SearchOptions = optimset('MaxFunEvals',10^6);  % By the way, SearchOptions can be used
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals,SearchOptions);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);
disp('The constraint implies that the fitted distributions should have equal means and std devs:');
fprintf('The estimated mean and standard deviation of Dist1 are %f and %f.\n',Dists{1}.Mean,Dists{1}.SD);
fprintf('The estimated mean and standard deviation of Dist2 are %f and %f.\n',Dists{2}.Mean,Dists{2}.SD);

%% Example: Constrain one distribution to be a mixture of two others.
% In this example we fit 2 normal distributions separately, but we also
% have a third distribution that is a mixture of the first two.
% 
% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

%%
% First generate a fake dataset (see **NOTE**).
n1=200; n2=190; n3=500;  % arbitrary sample sizes
TrueDist1 = Normal(100,10);
TrueDist2 = Normal(150,10);
TrueDist3 = Mixture(0.35,TrueDist1,.65,TrueDist2);
Datasets{1} = TrueDist1.Random(n1,1);   % Data from distributions.
Datasets{2} = TrueDist2.Random(n2,1);
Datasets{3} = TrueDist3.Random(n3,1);
%%
% Note that the constraint is truly satisfied because distribution is a mixture
% of distributions 1 and 2.
% Now we are done generating the fake datasets.
%
% Set some parameters defining the desired fit:
Dists{1} = Normal;
Dists{2} = Normal;
%%
% The parameters in the next line are not actually used but they must
% be specified so that the mixture distribution can be initialized properly.
Dists{3} = Mixture(.5,Normal(100,10),Normal(100,10));
sErrorFn = '-LnLikelihood';   % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [.5,100,10,150,15];      % Starting values for the parameters expected by DemoFitConstraintFn6.
ConstraintFn = @DemoFitConstrainedFns.Fn6;

%%
% Now perform the actual fit:
SearchOptions = optimset('MaxFunEvals',10^6);  % By the way, SearchOptions can be used
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals,SearchOptions);
fprintf('The fitted distributions are %s, %s,\n  and %s.',Dists{1}.StringName,Dists{2}.StringName,Dists{3}.StringName);

