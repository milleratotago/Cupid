% Demo of clOutlierModel
% Running the whole script takes ~10 minutes

prOutlier = 0.1;   % Probability that a true score is affected by the outlier process.

%% Example 1:  A single true distribution, sometimes shifted by an outlier process.

TrueDist = Normal(400,40);       % The distribution of true scores
ConType = clOutlierModel.Shift;  % In this model, contamination adds a random value to the true score
Contam = Normal(400,10);         % This is the distribution of the amount added

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.ObsDists{1} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(200,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)

%% Example 2:  Two different true distributions, each sometimes shifted by an outlier process.

TrueDist = {Normal(400,40), Normal(450,45)};       % CELL ARRAY with he distributions of true scores in two conditions
ConType = clOutlierModel.Shift;  % In this model, contamination adds a random value to the true score
Contam = Normal(400,35);         % This is the distribution of the amount added

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.TrueDist{2}, ...
    ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(500,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)


%% Example 3: 2 conditions, contamination by stretching.
% (all RVs must be strictly positive because CUPID's Product RV does not all 0 or negative values):

prOutlier = 0.2;   % Probability that a true score is affected by the outlier process.

TrueDist = {Beta(10,20), Beta(12,22)};
Contam = TriangularCW(1.4,0.2);
ConType = clOutlierModel.Stretch;
ModelWithContam = clOutlierModel(TrueDist,prOutlier,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.TrueDist{2}, ...
    ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(500,1);

% Example of parameter search (this takes 5-10 minutes, because Product is slow):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)


%% Example 4: 2 conditions, contamination by replacement:

prOutlier = 0.15;   % Probability that a true score is affected by the outlier process.

TrueDist = {RNGammaMS(400,40), RNGammaMS(450,45)};
Contam = RNGammaMS(1000,100);
ConType = clOutlierModel.Replace;

ModelWithContam = clOutlierModel(TrueDist,prOutlier,Contam,ConType);
X = ModelWithContam.Random(2000,1);

% Fix the means of the true dist, just for an illustration of fixing parameters.
parmCodes = ModelWithContam.DefaultParmCodes;
parmCodes([1 3]) = 'f';

% This is very fast, both because some parameters are fixed and because
% replacement is a simple mixture distribution without the Convolution
% or Product components used with Shift & Stretch.
[EndingVals,fval,exitflag,output] = EstML(ModelWithContam,X,parmCodes)

%% Example 5: 2 conditions, contamination by replacement, iterative estimation
%{

EstMLiter abandoned
prOutlier = 0.05;   % Probability that a true score is affected by the outlier process.

TrueDist = {RNGammaMS(400,40), RNGammaMS(450,45)};
Contam = RNGammaMS(1000,100);
ConType = clOutlierModel.Replace;

ModelWithContam = clOutlierModel(TrueDist,prOutlier,Contam,ConType);
X = ModelWithContam.Random(2000,1);

% Change the model to see if it can recover the right prOutlier!
ModelWithContam = clOutlierModel(TrueDist,2*prOutlier,Contam,ConType);

% This is very fast, both because some parameters are fixed and because
% replacement is a simple mixture distribution without the Convolution
% or Product components used with Shift & Stretch.
[EndingVals,fval,exitflag,output] = EstMLiter(ModelWithContam,X)
%}

%% Example 6: Fit with a known, fixed true distribution,
% only adjusting outlier probability and outlier distribution.

% Generate random CDFs with extra high-end outliers:
DataDist = Mixture(0.97,Uniform(0,1),Uniform(.99,1));
Observed = cell(2,1);
for iCond=1:2
    Observed{iCond} = DataDist.Random(10000,1);
end

% Model has Uniform(0,1) in all conditions (CDFs)
StartprOutlier = 0.05;   % Probability that a true score is affected by the outlier process.
TrueDist = {Uniform(0,1), Uniform(0,1)};
Contam = UniformCW(0.995,0.01);  % high-end outliers; use CW to avoid min>max in parameter searching
ConType = clOutlierModel.Replace;

ModelWithContam = clOutlierModel(TrueDist,StartprOutlier,Contam,ConType);

% fix true dist parms; free prOutlier & Contam parms:
ParmCodes = 'ffffrrr';
[EndingVals,fval,exitflag,output] = EstML(ModelWithContam,Observed,ParmCodes)

%% Example: 
% 2 conds: Identical true & contam dists but different probabilities.

TrueDist = Normal(400,40);       % The distribution of true scores
ConType = clOutlierModel.Replace;  % In this model, contamination replaces values of the true score
Contam = Normal(500,10);         % This is the distribution of the replacement value

prOutlier2 = [0.10, 0.20];
% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier2,Contam,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.ObsDists{1} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(2000,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X)


%% clOutlierModel 2 True Dists, Shift

TrueDist = {Uniform(0,1), Uniform(0,0.9)};  % The distributions of true scores
ConType = clOutlierModel.Shift;             % In this model, contamination add to values of the true score
ContamDists = UniformCW(1.9,0.009);         % This is the distribution of the additive values

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier,ContamDists,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
ModelWithContamb = clOutlierModel(TrueDist,0.22,ContamDists,ConType);
X = ModelWithContamb.Random(5000,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X,'ffffrrr')

figure;
for iCond=1:2
    subplot(1,2,iCond);
    plotObsPred(X{iCond},ModelWithContam.ObsDists{iCond});
end


%% clOutlierModel 2 True dists, Replace, 2 Contam dists

TrueDist = {Uniform(0,1), Uniform(0,0.9)};       %#ok<*UNRCH> % The distribution of CDFs of true scores
ConType = clOutlierModel.Replace;  % In this model, contamination replaces values of the true score
ContamDists = {UniformCW(0.9,0.2), UniformCW(0.95,0.1)};         % This is the distribution of the replacement value

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier,ContamDists,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(2000,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X,'ffffrrrrr')

%%  2 conditions, one true dist, two contam dists, replace

TrueDist = Uniform(0,1);       % The distribution of CDFs of true scores
ConType = clOutlierModel.Replace;  % In this model, contamination replaces values of the true score
ContamDists = {UniformCW(0.9,0.2), UniformCW(0.95,0.1)};         % This is the distribution of the replacement value

% Instantiate the model for the possibly-contaminated scores:
ModelWithContam = clOutlierModel(TrueDist,prOutlier,ContamDists,ConType);

% View the PDFs of the true scores and the observed (including contamination) scores:
figure; plotDists({ ModelWithContam.TrueDist{1}, ModelWithContam.ObsDists{1}, ModelWithContam.ObsDists{2} },1);

% Generate some random observed scores.
% These are a mixture of True and Contaminated (shifted true) scores
X = ModelWithContam.Random(200,1);

% Example of parameter search (this takes a minute or two):
[EndingVals,fval,exitflag,output] = ModelWithContam.EstML(X,'ffrrrrr')


