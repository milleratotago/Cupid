% messing about with exgamma & exwaldmsm

clear all
close all
clc

k = 5;
rate = 0.01;
tau = 100;
g = RNGamma(k,rate);
[g.Mean, g.Variance]
e = Exponential(1/tau);
[e.Mean, e.Variance]

eg = ExGamma(k,rate,1/tau);
[eg.Mean, eg.Variance]
x = eg.Random(5000,1);

a = eg.StartParmsMLE(x)
eg.ResetParms(a);
eg.PlotDens;
hold on
histogram(x,'normalization','pdf')

ew = ExWaldMSM(k/rate,sqrt(k/rate^2),tau);
xw = ew.Random(5000,1);
a = ew.StartParmsMLE(xw)
ew.ResetParms(a);
figure
ew.PlotDens;
hold on
histogram(xw,'normalization','pdf')
