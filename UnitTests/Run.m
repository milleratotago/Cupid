function Run(ThisCase)

%% Main script to run a set of tests

% Defaults:
% clc
% clear all
% close all   % old figs
import matlab.unittest.TestSuite  % avoids having to specify this as a qualifier on its method names:
runner = matlab.unittest.TestRunner.withTextOutput();
global WantPlots
global GlobalSkipAllEst
DiaryRoot = '00';

Continuous = false;
Discrete = false;
Derived = false;
Fast = false;
Slow = false;


% **************** OPTIONS START HERE ****************
% Note: Further speed options are available in utGeneric.m
% at the points marked *** CONTROL SPEED HERE ***

%  testCase.Dist.SearchOptions.Display = 'final';

% WantPlots = true;
WantPlots = false;

% GlobalSkipAllEst = true;
GlobalSkipAllEst = false;
if GlobalSkipAllEst
    DiaryRoot = [DiaryRoot 'GSE'];
end

Parallel = false;
% Parallel = true;

% WantBeep = true;
WantBeep = false;

warning off backtrace;   % Display 1-line warning messages.
% warning on backtrace;    % Display full warning messages.

switch ThisCase
    case 1
        Continuous = true;        Fast = true;
    case 2
        Continuous = true;        Slow = true;
    case 3
        Derived = true;          Fast = true;
    case 4
        Derived = true;          Slow = true;
    case 5
        Discrete = true;         Fast = true;
    case 6
        Discrete = true;         Slow = true;
end


% **************** OPTIONS END HERE ****************

tests = [];
results = [];

time0 = tic;

if Continuous && Fast
    StartDiary([DiaryRoot 'CntFas']);
    [t0, r0] = ClassTest(?utBeta,               Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utCauchy,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utChi,                Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utChiSq,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utCosine,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utDblMon,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExGauss,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExGauMn,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExGauRatio,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExponential,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExponenMn,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExpSum,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExtrVal1,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExtrVal2,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExtrValGen,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExWald,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExWaldMn,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExWaldMSM,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utF,                  Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utGamma,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utGenNor1,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utGenNor2,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utHyperbolicTan,      Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utJohnsonSB,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utJohnsonSU,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLaplace,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLogistic,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLognormal,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLognormalMS,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utNakaRush,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utNormal,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utPareto,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utQuantal,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utQuick,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRayleigh,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTriangular,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRNGamma,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRNGammaMn,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRNGammaMS,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRosin,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utStudRng,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utt,                  Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTriangularCW,       Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utUniform,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utUniformCW,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utUniGap,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utWald2,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utWeibull,            Parallel,runner); tests = [tests t0]; results = [results r0];
    diary('off');
end

if Continuous && Slow
    StartDiary([DiaryRoot 'CntSlo']);
    [t0, r0] = ClassTest(?utChiSqNoncentral,    Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExpSumT,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utFNoncentral,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utr,                  Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utrNoncentral,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRecinormal,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utSkewNor,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?uttNoncentral,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?uttPowerEst,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTriangularG,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTriangularGCWP,     Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utVonMises,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utWald,               Parallel,runner); tests = [tests t0]; results = [results r0];
    diary('off');
end

if Discrete && Fast
    StartDiary([DiaryRoot 'DisFas']);
    [t0, r0] = ClassTest(?utBinomial,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utBinomialMixe,       Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utGeometric,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utPoisson,            Parallel,runner); tests = [tests t0]; results = [results r0];
    % [t0, r0] = ClassTest(?utRankDist,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utUniformInt,         Parallel,runner); tests = [tests t0]; results = [results r0];
    diary('off');
end

if Discrete && Slow
    StartDiary([DiaryRoot 'DisSlo']);
%    [t0, r0] = ClassTest(?utBinomial,           Parallel,runner); tests = [tests t0]; results = [results r0];
%    [t0, r0] = ClassTest(?utPoisson,            Parallel,runner); tests = [tests t0]; results = [results r0];
    diary('off');
end

if Derived && Fast
    % This section takes ~14 min non-parallel
    StartDiary([DiaryRoot 'DrvFas']);
    [t0, r0] = ClassTest(?utAddTrans,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utArcsinTrans,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utExpTrans,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utInverseTrans,       Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLinearTrans,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLogTrans,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLogitTrans,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLogisticTrans,      Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utLogitLogistic,      Parallel,runner); tests = [tests t0]; results = [results r0]; % Double trans
    [t0, r0] = ClassTest(?utMixture,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utMonotoneTrans,      Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utMultTrans,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utPhiTrans,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utPhiInvTrans,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utPowerTrans,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utProduct,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utRatio,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utSqrTrans,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utSqrtTrans,          Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTruncatedP,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTruncatedX,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTruncatedXlow,      Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utTruncatedXhi,       Parallel,runner); tests = [tests t0]; results = [results r0];
    diary('off');
end

if Derived && Slow
    % This section takes 330 min non-parallel
    StartDiary([DiaryRoot 'DrvSlo']);
    [t0, r0] = ClassTest(?utAttainP,            Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utConvolution,        Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utDifference,         Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utInfMix,             Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utMinBound,           Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utOrder,              Parallel,runner); tests = [tests t0]; results = [results r0];
    [t0, r0] = ClassTest(?utOrderIID,           Parallel,runner); tests = [tests t0]; results = [results r0];
    diary('off');
end

% Note: When running parallel, MATLAB reports xxx seconds testing time, but it appears to sum the
% times for different processors working in parallel.

s=evalc('disp(results)');
wantloc=strfind(s,'Totals:');
fprintf(['Summary of test result totals:\n' s(wantloc+8:end)]);
totalmins_elapsed = toc(time0) / 60


if (totalmins_elapsed > 5) && WantBeep
    DoneBeep;
end

return

%% This section rechecks any that failed

close all   % old figs

failedTests = tests([results.Failed]);
result2 = run(failedTests)



%%
Parallel = false;

if Parallel
    runner = TestRunner.withTextOutput();
    result2 = runInParallel(runner,failedTests)
else
    result2 = run(failedTests)
end



%% Extra examples of specifying particular named tests to be run

% To specify a particular test, specify the test name as an argument of run. E.g.:
% utn = utNormal;
% utn.run('Moments');

% runtests('utNormal','ParameterName','zero','ParameterName','thousand')
% runtests('utAddTrans','ParameterName','parm3')


end


function StartDiary(DiaryName)
DiaryName = [DiaryName '.err'];
if exist(DiaryName,'file')>0
    delete(DiaryName)  % Delete old version of err file
end
diary(DiaryName);
end
