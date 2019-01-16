%% Multivariate random number generation demo
% This demo illustrates how the class "RandGen" can be used to generate multivariate random numbers.
% The marginals of the different RVs can be any distributions known to CUPID,
% and the marginal distributions can all be different.
%
% In the first several examples, RandGen functions are used directly in a non-object-oriented (OO) fashion.
% The last examples show how to use RandGen in an OO fashion.
%
% IMPORTANT: RandGen uses some routines from ExtractNameVal at https://github.com/milleratotago/ExtractNameVal
% so you need to have that in addition to Cupid.

%% A simple non-OO example

NCases = 25;  % Number of random cases to be generated.

%%
% Specify the number of random variables (RVs)
% and their marginal distributions.
NRVs = 3;
RVs = cell(NRVs,1);
RVs{1} = Normal(0,1);   % Each RV can be any legal Cupid distribution.
RVs{2} = Exponential(.01);
RVs{3} = Triangular(20,50);

%%
% Specify the upper triangular of a matrix determining
% the correlations among the RVs.
RhoControllers = ...
    [ 1  0.4  0.2; ...
    0   1   0.6; ...
    0   0    1   ];
RhoControllers = RhoControllers + triu(RhoControllers,1)';  % Copy the upper triangular portion into the lower triangular.
%%
% IMPORTANT NOTE:  The numbers in this RhoController matrix
% are not the _actual_ correlations among the RVs.
% They are only monotonically related to those correlations.
% There is more information about this in a later example
% describing how to adjust the RhoControllers to get
% a desired true correlation matrix.

randoms = RandGen.GenRands(RVs,RhoControllers,NCases)

%%
% Note that NCases is a vector, which is somewhat different than the usual
% convention when generating random numbers.  For example,
u = rand(10);
%%
% produces a 10x10 matrix of random numbers.

%% Relation of correlations to RhoController values (non-OO).
% Continuing on the simple example, this section shows that the true correlations
% among these RVs do NOT match the numbers in RhoControllers.  For the example, this
% is done by generating a large number of cases so that the true correlations can be
% estimated very accurately.  That is, the number of cases is large enough
% that the mismatches between computed correlations and RhoController values
% are not just due to random error.
NCases = 5000;
randoms = RandGen.GenRands(RVs,RhoControllers,NCases);
obscorrs = corr(randoms)
%%
% Note that the observed correlations do not match the values in RhoControllers.

%% Adjusting RhoController values to get desired correlations (non-OO).
% So how do you get the specific target correlations that you want?
% For convenience, GenRands can be told to adjust the RhoControllers
% automatically to give you the target correlations that you want.
% This allows you to skip the FindRhoControllerMatrix step from the previous section.
TargetCorrs = ...
    [ 1  0.5  0.3; ...
    0   1   0.5; ...
    0   0    1   ]

%%
% The optional argument 'Adjust' tells GenRands to find the required RhoControllers itself.
% With that argument, you can also specify a value for NStepsApprox (the default is 200).
NCases = 250000;   % Generate a lot of random numbers to get good estimates of the correlations.
randoms = RandGen.GenRands(RVs,TargetCorrs,NCases,'Adjust','NSteps',500);

corr(randoms)  % Check that the correlations are in fact quite close to our targets.

%% An example with discrete random variables (non-OO)
% Nothing is really different except that the random variables are discrete.
NRVs = 2;
RVs = cell(NRVs,1);
RVs{1}=Poisson(23);   % The Poisson parameter is the mean
RVs{2}=Poisson(12);
TargetCorrs = [1  0.3; 0.3 1];  % Replace 0.3 with whatever correlation you want.
NCases = 10000;
r = RandGen.GenRands(RVs,TargetCorrs,NCases,'Adjust','NSteps',500);
corr(r)
figure; histogram(r(:,1));


%% Example of object-oriented use
% Define the random variables pretty much as before.
NRVs = 3;
RVs = cell(NRVs,1);
RVs{1} = RNGamma(4,.05);   % Each RV can be any legal Cupid distribution.
RVs{2} = Exponential(.01);
RVs{3} = Triangular(20,50);

TargetCorrs = ...
    [ 1 -0.2  0.8; ...
    0   1  -0.4; ...
    0   0    1   ]

%%
% Now create a new random number generator object.
% It needs the same information as before to specify the multivariate distribution.
mymultivar = RandGen(RVs,TargetCorrs,'Adjust','NSteps',500,'Histograms','Scattergrams');

disp('Here are the values of RhoControllers, in case you want to save them for later re-use:');
temp1 = mymultivar.RhoControllers
disp('FYI, here are the maximum & minimal attainable correlations with these RVs:');
tempMax = mymultivar.RhoMax
tempMin = mymultivar.RhoMin

%%
% Now call the object's function to generate some random numbers from this distribution.
% By default, this function will also produce a histogram for each variable
% and a scatterplot for each pair of variables.
r=mymultivar.Rands(10000);

%%
% Here is how you would turn off automatic plotting of the two figures.
mymultivar.WantHistograms = false;
mymultivar.WantScattergrams = false;
r2=mymultivar.Rands(10000);  % No figures are plotted.

%%
% When would you want to use the OO approach?  Its main advantage
% is that it is faster if you generate many sets of random numbers
% from the same multivariate distribution.  This is because
% the object "remembers" the adjusted RhoControllers values
% so that you can get repeated samples from the same multivariate
% distribution without going through the (slow) RhoController
% adjustment step each time.  Here is an example:

NSamples = 3;
NCases = 10;

for iSample = 1:NSamples
    ThisSample = mymultivar.Rands(NCases)
    % Some processing of ThisSample
end

%% An OO example with discrete RVs

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

%%
% Here are the values of RhoControllers, in case you want to save them for later re-use.
% mymultivar holds onto them, but you might want to save them in a file to speed up later runs.
temp1 = mymultivar.RhoControllers
disp('FYI, here are the maximum & minimal attainable correlations with these RVs:');
tempMax = mymultivar.RhoMax
tempMin = mymultivar.RhoMin

r=mymultivar.Rands(10000);

mymultivar.WantHistograms = false;
mymultivar.WantScattergrams = false;

%%
% Now verify that the random numbers from this distribution have the desired correlations.
NCases = 250000;
ThisSample = mymultivar.Rands(NCases);
corr(ThisSample)
for iRV = 1:NRVs
    figure;
    histogram(ThisSample(:,iRV));
end





