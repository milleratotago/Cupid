% ********** Commands from cupid.tex

%% Distributions associated with hypothesis testing procedures:
df=19;
pcurve=AttainP(tNoncentral(df,.4),t(df));
pcurve.PlotDens

%% User Interface
mydist1=Uniform(0,100);
mydist2=Mixture(0.65,Normal(0,1),0.35,Normal(1,1));
mydist2.PlotDens;
my1sd=mydist1.SD

% ************ Command from cupdists.tex

%% **************** Standard Distributions

mydist=Beta(2.3,3.2);                       mydist.PlotDens;
mydist=Cauchy(2.3,3.2);                     mydist.PlotDens;
mydist=ChiSq(20);                           mydist.PlotDens;
mydist=ChiSqNoncentral(20,0.4);             mydist.PlotDens;
mydist=Chi(20);                             mydist.PlotDens;
mydist=Cosine(100,20);                      mydist.PlotDens;
mydist=DblMon(100,4,20);                    mydist.PlotDens;
mydist=ExGauss(300,30,.01);                 mydist.PlotDens;
mydist=ExGauMn(300,30,100);                 mydist.PlotDens;
mydist=ExGauRatio(300,30,1/3);              mydist.PlotDens;
mydist=Exponential(.05);                    mydist.PlotDens;
mydist=ExponenMn(20);                       mydist.PlotDens;
mydist=ExpSum(.01,.02);                     mydist.PlotDens;
mydist=ExpSumT(.01,.02,300);                mydist.PlotDens;
% mydist=ExpoNo(100,20);                      mydist.PlotDens;
mydist=ExtrVal1(100,10);                    mydist.PlotDens;
mydist=ExtrVal2(100,20,23);                 mydist.PlotDens;
mydist=ExtrValGen(100,20,.5);               mydist.PlotDens;
mydist=ExWald(1,2,20,.3);                   mydist.PlotDens;
mydist=ExWaldMn(1,2,20,30);                 mydist.PlotDens;
mydist=ExWaldMSM(100,20,30);                mydist.PlotDens;
mydist=F(4,20);                             mydist.PlotDens;
mydist=FNoncentral(4,20,3);                mydist.PlotDens;
mydist=Gamma(5,.2);                         mydist.PlotDens;
% mydist=Geary(50);                           mydist.PlotDens;
% mydist=GenErr(100,20,.5);                   mydist.PlotDens;
mydist=HyperbolicTan(.3);                   mydist.PlotDens;
mydist=JohnsonSB(100,20,.4,1.4);            mydist.PlotDens;
mydist=JohnsonSU(100,20,.4,1.4);            mydist.PlotDens;
% mydist=Kolmog(100);                         mydist.PlotDens;
mydist=Laplace(100,20);                     mydist.PlotDens;
% mydist=Lilliefors(100);                     mydist.PlotDens;
mydist=Logistic(100,20);                    mydist.PlotDens;
mydist=Lognormal(2.3,.1);                   mydist.PlotDens;
mydist=LognormalMS(20,5);                   mydist.PlotDens;
mydist=NakaRush(10);                        mydist.PlotDens;
mydist=Normal(100,10);                      mydist.PlotDens;
% mydist=NormalLowerHalf(100,10);             mydist.PlotDens;
% mydist=NormalUpperHalf(100,10);             mydist.PlotDens;
mydist=Pareto(25,8);                        mydist.PlotDens;
mydist=Quantal(16);                         mydist.PlotDens;
mydist=Quick(8,12);                         mydist.PlotDens;
mydist=r(20);                               mydist.PlotDens;
mydist=rNoncentral(20,.5);                  mydist.PlotDens;
mydist=Rayleigh(12);                        mydist.PlotDens;
mydist=Recinormal(.005,.001);               mydist.PlotDens;
mydist=RNGamma(3.2,.025);                   mydist.PlotDens;
mydist=RNGammaMn(3.2,200);                  mydist.PlotDens;
mydist=RNGammaMS(302,120);                  mydist.PlotDens;
mydist=Rosin(10,1.5);                       mydist.PlotDens;
mydist=SkewNor(-20,3,-5);                   mydist.PlotDens;
% mydist=StudRng(20,4);                       mydist.PlotDens;
mydist=t(20);                               mydist.PlotDens;
mydist=tNoncentral(20,.4);                  mydist.PlotDens;
mydist=Triangular(0,1);                     mydist.PlotDens;
mydist=TriangularCW(0,1);                   mydist.PlotDens;
mydist=TriangularG(0,.8,1);                 mydist.PlotDens;
mydist=TriangularGCWP(0,100,.8);            mydist.PlotDens;
mydist=Uniform(0,100);                      mydist.PlotDens;
mydist=UniformCW(50,50);                    mydist.PlotDens;
% mydist=UniGap(25);                          mydist.PlotDens;
mydist=VonMises(.85,.1);                    mydist.PlotDens;
mydist=Wald(.20,2,50);                      mydist.PlotDens;
mydist=Wald2(.20,50);                       mydist.PlotDens;
mydist=Weibull(10,2.5,0);                   mydist.PlotDens;

%% Transformation distributions
mydist=AddTrans(Uniform(.5,1),10);          mydist.PlotDens;
mydist=ArcsinTrans(Uniform(.5,1));          mydist.PlotDens;
mydist=ExpTrans(Uniform(.5,1));             mydist.PlotDens;
mydist=InverseTrans(Uniform(.5,1));         mydist.PlotDens;
mydist=LinearTrans(Uniform(.5,1),2,10);     mydist.PlotDens;
mydist=LogTrans(Uniform(0.5,1));            mydist.PlotDens;
mydist=MultTrans(Uniform(.5,1),2);          mydist.PlotDens;
mydist=PhiTrans(Uniform(-1,1));             mydist.PlotDens;
mydist=PhiInvTrans(Beta(3,3));              mydist.PlotDens;
mydist=PowerTrans(Uniform(.5,1),2);         mydist.PlotDens;
mydist=SqrTrans(Uniform(.5,1));             mydist.PlotDens;
mydist=SqrtTrans(Uniform(.5,1));            mydist.PlotDens;


% **************** Derived distributions

%% Convolution distributions

mydist=Convolution(Normal(0,1),Uniform(0,1));
mydist.PlotDens;

mydist=Convolution(Normal(100,50),Convolution(Uniform(0,100),Gamma(3,0.01)));
% This distribution takes too long to plot, so we will just get a few
% easy-to-compute values.
% mydist.PlotDens;
mydist.PDF(1)
mydist.Mean
mydist.SD


%% Difference distributions

mydist=Difference(Uniform(0,1),Uniform(0,1));
mydist.PlotDens;


%% Mixtures and Infinite Mixtures

mydist=Mixture(0.5,Normal(0,1),0.5,Uniform(0,1));
mydist.PlotDens;

mydist=InfMix(Normal(0,5),Uniform(10,20),1);
mydist.PlotDens;

mydist=InfMix(Normal(0,5),Uniform(10,20),2);
mydist.PlotDens;

mydist=InfMix(InfMix(Normal(0,5),Uniform(0,2),1),Uniform(4,6),2);
mydist.PlotDens;


%% Product and Ratio

mydist=Product(Uniform(0,1),Uniform(0,1));
mydist.PlotDens;

mydist=Ratio(Uniform(0,1),Uniform(1,2));
mydist.PlotDens;


%% Truncated distributions

mydist=TruncatedX(Normal(0,1),-1,1);
mydist.PlotDens;

mydist=TruncatedX(BasisDistribution(Parms),Min,Max);
mydist.PlotDens;

mydist=TruncatedP(BasisDistribution(Parms),0.05,0.95);
mydist.PlotDens;


%% Order statistics

mydist=Order(2,Normal(0,1),Uniform(0,1),Exponential(1));
mydist.PlotDens;

mydist = OrderIID(3,10,Exponential(.1));
mydist.PlotDens;


%% MinBound distributions

xvalues = 0.01:.01:.99;
mydist1 = Beta(5.4,10.5);
mydist2 = Uniform(0,1);
mydistmin = MinBound(mydist1,mydist2);
mydist1cdf = mydist1.CDF(xvalues);
mydist2cdf = mydist2.CDF(xvalues);
mydistmincdf = mydistmin.CDF(xvalues);
figure;
plot(xvalues,mydistmincdf);
hold on;
plot(xvalues,mydist1cdf);
plot(xvalues,mydist2cdf);
legend({'MinBound','Beta','Uniform'});


% **************** Distributions associated with hypothesis testing

%% AttainP

mydist=AttainP(Normal(1,1),Normal(0,1));
mydist.PlotDens;

%% tPowerEst

mydist=tPowerEst(.5,1,.05,20);
mydist.PlotDens;


%% Examples of recursive distribution definition

mydist=TruncatedX(Mixture(.5,Normal(0,1),.5,OrderIID(4,5,Normal(0,1))),-1,1);
mydist.PlotDens;

mydist=Convolution(OrderIID(1,100,Normal(0,1)),Normal(1,1));
mydist.PlotDens;


% **************** Commands describing functions

%% Distribution functions

x=0.33;
mydist.PDF(x)
mydist.CDF(x)
mydist.OMCDF(x)
mydist.TwoTailProb(x)
mydist.Hazard(x)

%% Moment-based functions
mydist.Mean
mydist.SD
mydist.Variance
mydist.CV
mydist.Skewness
mydist.RawSkewness
mydist.Kurtosis

%% Percentile functions
mydist.InverseCDF(P)
mydist.Median
mydist.SIQR
mydist.PctileSkew(P)
mydist.MMMSD

%% Integrals
mydist.RawMoment(Power)
mydist.ConditionalRawMoment(Min,Max,Power)
mydist.IntegralXtoNxPDF(Min,Max,Power)
% mydist.IntegrateOverP=true;
mydist.IntegralX_CToNxPDF(Min,Max,C,Power)
mydist.CDFIntegral(Min,Max,Power)
mydist.MGF(Theta)

%% Goodness-of-fit

mydist = Beta(1.4,1.5);
B = .1:.1:1;
O = [.1 .11 .09 .1 .12 .08 .1 .1 .09 .11];
GoF = mydist.GofFChiSq(B, O)

X = mydist.Random(10);
mydist.LnLikelihood(X)

%% Random numbers
mydist.Random
mydist.Random(6)
mydist.RndLnLikelihood(6)



% **************** Other Functionality

%% MakeBinSet
a=Uniform(0,1).MakeBinSet(10,true)
Normal(0,1).MakeBinSet(5,false)

%% FnAfterReset
Uniform(0,1).FnAfterReset(2,1:10,'Mean')

%% PrXGTY
thisp = PrXGTY(Normal(0,1),Uniform(0,1))

% **************** Estimation

%% Estimate From Moments

mydist=RNGamma(3,.4);
mydist.EstMom([100,40000])


%% Estimate From ChiSquare

mydist = Beta(1.4,1.5);
B = .1:.1:1;
O = [.1 .11 .09 .1 .12 .08 .1 .1 .09 .11];
mydist.EstChiSq(B,O)


%% Estimate From Percentiles

mydist = Beta(2.4,1.5);
B = .1:.1:.9;  % Values of the beta RV
P = [.1 .21 .29 .30 .42 .48 .5 .7 .99];
mydist.EstPctile(B,P)


%% Maximum Likelihood Estimation

mydist = Lognormal(2.4,1.5);
X = mydist.Random(100,1)';  % Generate 100 random numbers
mydist.EstML(X)


%% Interval Probability Estimation

mydist = Lognormal(2.4,1.5);
mydist.EstPctBounds(2,10,.3)

%% ParmCodes

mydist=Beta(0.5,0.5);
mydist.EstMom([0.3,0.1],'rr')

mydist=Beta(0.5,0.5);
mydist.EstMom([0.3,0.1],'rf')

mydist=RNGamma(3.5,10.5);
mydist.EstML([0.34162, 0.52264, 0.35699, 0.40554, 0.34145, 0.37642],'ii')

mydist=TruncatedX(Normal(0,1),-1,1);
mydist.EstMom([0.3,0.22],'rrff')

mydist=Mixture(0.6,Normal(0,1),0.4,Exponential(1));
mydist.EstML([4,1.1,3.2,-0.3,0.8],'rfff')


%% A Shortcut for Changing Parameter Values

mydist=Convolution(TruncatedX(Normal(0,1),-1,1),OrderIID(3,5,Exponential(.10)));
mydist.StringName

mydist.ResetSomeParms(7,.15);
mydist.StringName

mydist.ResetSomeParms(1,0.033,7,.15)
mydist.StringName


%% Spline Approximations

mydist = InfMix(Normal(0,1),RNGamma(3,.2),1);
tic
mydist.PlotDens
toc

mydist = InfMix(Normal(0,1),RNGamma(3,.2),1);
tic
mydist.UseSplinePDFOn(50);
mydist.PlotDens
toc

mypart = Convolution(Uniform(0,100),Gamma(3,0.01));
mypart.UseSplinePDFOn(100);
mydoubleconv = Convolution(Normal(0,50),mypart);
mydoubleconv.PlotDens;


%% Programming with CUPID Handle Objects

MyNorm = Normal(0,1);
MyWinner = Order(1,MyNorm,Exponential(1));
[MyWinner.Mean MyWinner.SD]   % MyWinner has a mean & sd of -0.16052 and 0.8224.

MyNorm.ResetSomeParms(1,2);   % Change MyNorm's mean to 2
[MyWinner.Mean MyWinner.SD]   % MyWinner now has a mean & SD of 0.78103 and 0.70097


