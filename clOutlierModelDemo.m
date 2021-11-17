% Demo of clOutlierModel

prOutlier = 0.1;   % Probability that a true score is affected by the outlier process.

%% Example 1:  A single true distribution, sometimes shifted by an outlier process.

TrueDist = Normal(400,40);       % The distribution of true scores
ConType = clOutlierModel.Shift;  % In this model, contamination adds a random value to the true score
Contam = Normal(400,10);         % This is the distribution of the amount added

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDists{1}, ModelWithContam.ObsDists{1} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(200,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)

%% Example 2:  Two different true distributions, each sometimes shifted by an outlier process.

TrueDists = {Normal(400,40), Normal(450,45)};       % CELL ARRAY with he distributions of true scores in two conditions
ConType = clOutlierModel.Shift;  % In this model, contamination adds a random value to the true score
Contam = Normal(400,35);         % This is the distribution of the amount added

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDists,prOutlier,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDists{1}, ModelWithContam.TrueDists{2}, ...
    ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(500,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)


%% Example 3: 2 conditions, contamination by stretching.
% (all RVs must be strictly positive because CUPID's Product RV does not all 0 or negative values):

prOutlier = 0.2;   % Probability that a true score is affected by the outlier process.

TrueDists = {Beta(10,20), Beta(12,22)};
Contam = Triangular(1.3,1.5);
ConType = clOutlierModel.Stretch;
ModelWithContam = clOutlierModel(TrueDists,prOutlier,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDists{1}, ModelWithContam.TrueDists{2}, ...
    ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(500,1);

% Example of parameter search (this takes 5-10 minutes, because Product is slow):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)



%% Example 4: 2 conditions, contamination by replacement:

prOutlier = 0.15;   % Probability that a true score is affected by the outlier process.

TrueDists = {RNGammaMS(400,40), RNGammaMS(450,45)};
Contam = RNGammaMS(1000,100);
ConType = clOutlierModel.Replace;

ModelWithContam = clOutlierModel(TrueDists,prOutlier,Contam,ConType);
X = ModelWithContam.Random(2000,1);

% Fix the means of the true dist, just for an illustration of fixing parameters.
parmCodes = ModelWithContam.DefaultParmCodes;
parmCodes([1 3]) = 'f';

% This is very fast, both because some parameters are fixed and because
% replacement is a simple mixture distribution without the Convolution
% or Product components used with Shift & Stretch.
[EndingVals,fval,exitflag,output] = EstML(ModelWithContam,X,parmCodes)

%% Example 5: 2 conditions, contamination by replacement, iterative estimation

prOutlier = 0.05;   % Probability that a true score is affected by the outlier process.

TrueDists = {RNGammaMS(400,40), RNGammaMS(450,45)};
Contam = RNGammaMS(1000,100);
ConType = clOutlierModel.Replace;

ModelWithContam = clOutlierModel(TrueDists,prOutlier,Contam,ConType);
X = ModelWithContam.Random(2000,1);

% Change the model to see if it can recover the right prOutlier!
ModelWithContam = clOutlierModel(TrueDists,2*prOutlier,Contam,ConType);

% This is very fast, both because some parameters are fixed and because
% replacement is a simple mixture distribution without the Convolution
% or Product components used with Shift & Stretch.
[EndingVals,fval,exitflag,output] = EstMLiter(ModelWithContam,X)
