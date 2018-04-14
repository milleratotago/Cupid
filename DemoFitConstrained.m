% Demo script for fitting contrained distribution(s).

% **NOTE**: In each example, the first step is to generate a fake dataset
% so that there are data available for fitting. Presumably you would already
% have the data that you wanted to fit, so you would skip this step.

%% 
disp('Example 0: Maximum likelihood fit of a normal distribution with the constraint sigma>0.1*mu.')

% First generate a fake dataset (see **NOTE**).
n1=100;  % arbitrary sample size
Datasets{1} = normrnd(55,4.5,n1,1);
% Note that the constraint is not truly satisfied because the true sigma is less than 5.5.

% Set some parameters defining the desired fit:
Dists{1} = Normal(45,4);       % The distribution to be fit.  These parameter values are not actually used.
sErrorFn = '-LnLikelihood';    % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [40, sqrt(5)];  % Starting values for the two parameters accepted by DemoFitConstraintFn0.
ConstraintFn = @DemoFitConstraintFn0;   % This user-supplied function converts fminsearch's suggested parameter values
                                        % into suggested parameter values for the distributions that are being fit.
                                        % Look at the function to see what it does.
% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distribution is ' Dists{1}.StringName]);

% What would the fit have been without the constraint?
nocon = Dists{1}.EstML(Datasets{1});
disp(['In comparison, without the constraint the estimated distribution is:' nocon]);

disp(' ');

%% 
disp('Example 1: Maximum likelihood fit of two normal distributions with the constraint that they have the same sigma.');

% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

% First generate a fake dataset (see **NOTE**).
n1=100; n2=90;  % arbitrary sample sizes
Datasets{1} = normrnd(10,10,n1,1);
Datasets{2} = normrnd(20,12,n2,1);
% Note that the constraint is not truly satisfied because the true sigmas are not equal.

% Set some parameters defining the desired fit:
Dists{1} = Normal(15,11);  % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = Normal(15,11);
sErrorFn = '-LnLikelihood';  % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [15, 15, sqrt(11)];      % Starting values for the three parameters accepted by DemoFitConstraintFn1.
ConstraintFn = @DemoFitConstraintFn1;   % This user-supplied function converts fminsearch's suggested parameter values
                                        % into suggested parameter values for the distributions that are being fit.

% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);

disp(' ');

%% 
disp('Example 2: Maximum likelihood fit of three distributions: a normal, an exponential, and an ex-Gaussian, with the');
disp('constraint that the parameters of the ex-Gaussian are equal to those of the normal and exponential.');

% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

% First generate a fake dataset (see **NOTE**).
n1=100; n2=90; n3 = 130;  % arbitrary sample sizes
Datasets{1} = normrnd(300,30,n1,1);
Datasets{2} = exprnd(100,n2,1);
Datasets{3} = normrnd(320,32,n3,1) + exprnd(100,n3,1);
% Note that the constraint is not truly satisfied because the parameters of the exGaussian
% are not identical to those of the normal and exponential.

% Set some parameters defining the desired fit:
Dists{1} = Normal(15,11);       % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = ExponenMn(100);      % Use the exponential whose parameter is the mean, not rate.
Dists{3} = ExGauMn(15,11,100);  % Use the exGaussian whose exponential parameter is the mean, not rate.
sErrorFn = '-LnLikelihood';     % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [290, sqrt(23), sqrt(100)];  % Starting values for the three parameters accepted by DemoFitConstraintFn2.
ConstraintFn = @DemoFitConstraintFn2;       % This user-supplied function converts fminsearch's suggested parameter values
                                            % into suggested parameter values for the distributions that are being fit.

% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ', ' Dists{2}.StringName ', and ' Dists{3}.StringName]);

disp(' ');

%%
disp('Example 3: ChiSquare fit of two normal distributions with the constraint that they have the same sigma.');
disp('This example requires construction of more complicated Datasets, because the GofFChiSq error function');
disp('requires both bin specifications and counts as input data.');

% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

% First generate a fake dataset (see **NOTE**).
n1=400; n2=590;  % arbitrary sample sizes
Scores1 = normrnd(10,10,n1,1);
Scores2 = normrnd(20,12,n2,1);
% Note that the constraint is not truly satisfied because the true sigmas are not equal.
%
% Use histcounts to bin and count the scores, since GofFChiSq works
% with bins and probabilities, not raw scores.
[counts1,bins1] = histcounts(Scores1);
[counts2,bins2] = histcounts(Scores2);
%
% Make the datasets that will be passed to GofFChiSq.
% Each dataset specifies GofFChiSq's parameters of BinUpperBounds and BinProbs, in that order.
Datasets{1}{1} = bins1(2:end);
Datasets{1}{2} = counts1/sum(counts1);
Datasets{2}{1} = bins2(2:end);
Datasets{2}{2} = counts2/sum(counts2);
% Now we are done generating the fake datasets.

% Set some parameters defining the desired fit:
Dists{1} = Normal(15,11);  % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = Normal(15,11);
sErrorFn = 'GofFChiSq';  % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [15, 15, sqrt(11)];      % Starting values for the three parameters accepted by DemoFitConstraintFn1.
ConstraintFn = @DemoFitConstraintFn1;   % This user-supplied function converts fminsearch's suggested parameter values
                                        % into suggested parameter values for the distributions that are being fit.

% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);

disp(' ');


%%
disp('Example 4: Fit two logistic distributions to two 2AFC psychometric functions with');
disp('the constraint described by Bausenhart et al, Behav Res (2012) 44:1157–1174');
disp('DOI 10.3758/s13428-012-0207-z');
disp('This example also requires construction of complicated Datasets, because');
disp('the mAFCProbitLnLikelihood error function requires multiple inputs.');

% Clear these variables so that MATLAB will not be confused by their previous types.
clear Dists;
clear Datasets;

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
% Now we are done generating the fake datasets.

% Set some parameters defining the desired fit:
Dists{1} = Logistic(50,10);  % The distributions to be fit.  These parameter values are not actually used.
Dists{2} = Logistic(50,10);
sErrorFn = '-mAFCProbitLnLikelihood';   % The name of the error function to be minimized; see FitConstrained.m for a list of the options.
StartingVals = [15, 15, sqrt(11)];      % Starting values for the three parameters accepted by DemoFitConstraintFn4.
ConstraintFn = @DemoFitConstraintFn4;   % This user-supplied function converts fminsearch's suggested parameter values
                                        % into suggested parameter values for the distributions that are being fit.

% Now perform the actual fit:
[Dists, ErrScores] = FitConstrained(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals);
disp(['The fitted distributions are ' Dists{1}.StringName ' and ' Dists{2}.StringName]);
disp('The constraint implies that the CDFs of these distributions evaluated at s=50 should sum to 1.0:');
disp(['The actual sum is ' num2str(Dists{1}.CDF(s) + Dists{2}.CDF(s))]);

