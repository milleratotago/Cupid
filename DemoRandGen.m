% Random number generation demos.

% IMPORTANT: RandGen uses some routines from ExtractNameVal at https://github.com/milleratotago/ExtractNameVal
% so you need to have that in addition to Cupid.

%% A simple example:

NCases = 25;

% Specify the number of random variables (RVs)
% and their marginal distributions.
NRVs = 3;
RVs = cell(NRVs,1);
RVs{1} = Normal(0,1);   % Each RV can be any legal Cupid distribution.
RVs{2} = Exponential(.01);
RVs{3} = Triangular(20,50);

% Specify the upper triangular of a matrix determining
% the correlations among the RVs.
RhoControllers = ...
    [ 1  0.4  0.2; ...
    0   1   0.6; ...
    0   0    1   ];
RhoControllers = RhoControllers + triu(RhoControllers,1)';  % Copy the upper triangular portion into the lower triangular.
% IMPORTANT NOTE:  The numbers in this RhoController matrix
% are not the _actual_ correlations among the RVs.
% They are only monotonically related to those correlations.
% There is more information about this in a later example
% describing how to adjust the RhoControllers to get
% a desired true correlation matrix.

randoms = RandGen.GenRands(RVs,RhoControllers,NCases)

%% Continuing on the simple example, this section shows that the true correlations
% among these RVs do NOT match the numbers in RhoControllers.  For the example, this
% is done by generating a large number of cases so that the true correlations can be
% estimated very accurately.  That is, the number of cases is large enough
% that the mismatches between computed correlations and RhoController values
% are not just due to random error.
NCases = 5000;
randoms = RandGen.GenRands(RVs,RhoControllers,NCases);
obscorrs = corr(randoms)
% Note that the observed correlations do not match the values in RhoControllers.

%% So how do you get the specific target correlations that you want?
% For convenience, GenRands can be told to adjust the RhoControllers
% automatically to give you the target correlations that you want.
% This allows you to skip the FindRhoControllerMatrix step from the previous section.
TargetCorrs = ...
    [ 1  0.5  0.3; ...
    0   1   0.5; ...
    0   0    1   ]

% The optional argument 'Adjust' tells GenRands to find the required RhoControllers itself.
% With that argument, you can also specify a value for NStepsApprox (the default is 200).
NCases = 250000;   % Generate a lot of random numbers to get good estimates of the correlations.
randoms = RandGen.GenRands(RVs,TargetCorrs,NCases,'Adjust','NSteps',500);

corr(randoms)  % Check that the correlations are in fact quite close to our targets.


%% Example of object-based use:

NRVs = 3;
RVs = cell(NRVs,1);
RVs{1} = RNGamma(4,.05);   % Each RV can be any legal Cupid distribution.
RVs{2} = Exponential(.01);
RVs{3} = Triangular(20,50);

TargetCorrs = ...
    [ 1 -0.2  0.8; ...
    0   1  -0.4; ...
    0   0    1   ]

mymultivar = RandGen(RVs,TargetCorrs,'Adjust','NSteps',500,'Histograms','Scattergrams');

disp('Here are the values of RhoControllers, in case you want to save them for later re-use:');
temp1 = mymultivar.RhoControllers
disp('FYI, here are the maximum & minimal attainable correlations with these RVs:');
tempMax = mymultivar.RhoMax
tempMin = mymultivar.RhoMin

r=mymultivar.Rands(10000);

mymultivar.WantHistograms = false;
mymultivar.WantScattergrams = false;

NSamples = 3;
NCases = 10;

for iSample = 1:NSamples
    ThisSample = mymultivar.Rands(NCases)
end

%% Example with discrete RVs:

NRVs = 3;
RVs = cell(NRVs,1);
RVs{1} = Poisson(4);   % Each RV can be any legal Cupid distribution.
RVs{2} = UniformInt(1,10);
RVs{3} = Binomial(20,.5);

TargetCorrs = ...
    [ 1 -0.2  0.8; ...
    0   1  -0.4; ...
    0   0    1   ]

mymultivar = RandGen(RVs,TargetCorrs,'Adjust','NSteps',500,'Histograms','Scattergrams');

disp('Here are the values of RhoControllers, in case you want to save them for later re-use:');
temp1 = mymultivar.RhoControllers
disp('FYI, here are the maximum & minimal attainable correlations with these RVs:');
tempMax = mymultivar.RhoMax
tempMin = mymultivar.RhoMin

r=mymultivar.Rands(10000);

mymultivar.WantHistograms = false;
mymultivar.WantScattergrams = false;

% % Now verify that these random numbers have the desired correlations.
NCases = 250000;
ThisSample = mymultivar.Rands(NCases);
corr(ThisSample)
for iRV = 1:NRVs
    figure;
    histogram(ThisSample(:,iRV));
end


%% An example with discrete random variables:

NRVs = 2;
RVs = cell(NRVs,1);
RVs{1}=Poisson(23);   % The Poisson parameter is the mean
RVs{2}=Poisson(12);
TargetCorrs = [1  0.3; 0.3 1];  % Replace 0.3 with whatever correlation you want.
mymultivar = RandGen(RVs,TargetCorrs,'Adjust','NSteps',500);
r=mymultivar.Rands(10000);  % r is the array of random numbers with desired marginals & correlation
corr(r)
histogram(r(:,1));

%% %%%%%%%%%%%%%%%%% Old Versions from here:

% NwJeff: Spline examples.  Add 'SplineAll' option?


% %% This section shows how to adjust the values in RhoControllers to get
% % any desired true correlations among RVs (to a good approximation).
%
% % Suppose that we would actually like to have these correlations:
% TargetCorrs = ...
%   [ 1  0.5  0.3; ...
%     0   1   0.5; ...
%     0   0    1   ];
%
% % Determine what RhoControllers values will give the desired target correlations.
% % This is done by approximating each distribution with NStepsApprox equally-spaced
% % (in percentiles) steps. (For example, with 100 steps the distribution would be
% % approximated using its values at the percentiles .5, 1.5, 2.5, ... 99.5.)
% NStepsApprox = 200;
% RhoControllers = RandGen.FindRhoControllerMatrix(RVs,TargetCorrs,NStepsApprox)
%
% % Now verify that these RhoController values give the desired correlations.
% NCases = 250000;
% randoms = RandGen.GenRands(RVs,RhoControllers,NCases);  % Generate a lot of random numbers using those RhoController values.
% corr(randoms)  % Check that the correlations are in fact quite close to our targets.
