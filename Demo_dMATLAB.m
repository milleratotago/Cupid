% Examples of using native MATLAB probability distribution objects within Cupid.
% In this demo file, the m* objects are MATLAB probability distributions,
% and the c* objects are Cupid distributions.

%% Set up and use the Rician distribution as an example:

mRician = makedist('Rician','s',2,'sigma',3);
cRician = dMATLABc(mRician,'rr',[0 0],[+inf +inf]);  % s & sigma are positive reals.

cRician.PlotDens;
[cRician.Mean, cRician.Variance]

% Adjust Rician parameters to give desired mean and variance of 5.
warning('off','stats:ncx2inv:NotConverge');  % There are a lot of warnings that NCX2INV does not converge.
warning('off','stats:ncx2inv:LastStep');
cRician.EstMom([5 5])

[cRician.Mean, cRician.Variance]

%% Form a convolution of two identical Ricians

cConv = Convolution(cRician,cRician);
% Note that this convolution has two copies of the same distribution object so they
% will necessarily always have same parameter values (even if parameters are adjusted).
[cConv.Mean, cConv.Variance]

% Adjust the parameters so that the convolution has a mean and variance of 11
cConv.EstMom([11 11],'ffrr')  % To save time, tell Cupid not to adjust the parameters of the first distribution,
                              % since they are actually identical to the second--it is all one distribution.
                              % Unfortunately, 'rrff' does not also work, because of how Cupid sets parameters.

[cConv.Mean, cConv.Variance]

cConv.PlotDens

%% Look at the distribution of the minimum of two different Ricians:

% For simplicity, reset the parameters of the Rician
cRician.ResetParms([2 3]);

% Make two new distributions so that they can have different parameter values than the original ones.
mRician2 = makedist('Rician','s',2.5,'sigma',3);
cRician2 = dMATLABc(mRician2,'rr',[0 0],[+inf +inf]);  % s & sigma are positive reals.

cOrder = Order(1,cRician,cRician2);
cOrder.PlotDens;
cOrder.CDF([2 4])

% Adjust the parameters of the second Rician so that the minimum has the target CDF(2)=.25 and CDF(5)=.75:
cOrder.EstPctile([2 4],[.25 .75],'fffrr')  % fix the ith order parameter and the two parameter values of the first Rician
cOrder.CDF([2 4])

