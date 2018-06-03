function Run(CaseNs)
    % Function to run a set of tests indicated by a vector of Ns 1-8.
    
    global WantPlots
    global GlobalSkipAllEst

    % **************** OPTIONS START HERE ****************
    % Note: Further speed options are available in utGeneric.m
    % at the points marked *** CONTROL SPEED HERE ***
    
    WantPlots = false;
    GlobalSkipAllEst = false;
    % WantBeep = true;
    WantBeep = false;
    
    % **************** OPTIONS END HERE ****************

    MarkerRoot = '00Results/00';
    if GlobalSkipAllEst
        MarkerRoot = [MarkerRoot 'GSE'];
    end
        
    for iCase=CaseNs
        
        Continuous = false;
        Discrete = false;
        DerivedCont = false;
        DerivedDisc = false;
        Fast = false;
        Slow = false;
        
        switch iCase
            case 1; Continuous  = true; Fast = true;
            case 2; Discrete    = true; Fast = true;
            case 3; DerivedCont = true; Fast = true;
            case 4; DerivedDisc = true; Fast = true;
            case 5; Continuous  = true; Slow = true;
            case 6; Discrete    = true; Slow = true;  % UNUSED
            case 7; DerivedCont = true; Slow = true;
            case 8; DerivedDisc = true; Slow = true;
        end
        
        time0 = tic;
        
        if Continuous && Fast
            SingleTest('Beta');
            SingleTest('Cauchy');
            SingleTest('Chi');
            SingleTest('ChiSq');
            SingleTest('Cosine');
            SingleTest('DblMon');
            SingleTest('ExGauss');
            SingleTest('ExGauMn');
            SingleTest('ExGauRatio');
            SingleTest('Exponential');
            SingleTest('ExponenMn');
            SingleTest('ExpSum');
            SingleTest('ExtrVal1');
            SingleTest('ExtrVal2');
            SingleTest('ExtrValGen');
            SingleTest('ExWald');
            SingleTest('ExWaldMn');
            SingleTest('ExWaldMSM');
            SingleTest('F');
            SingleTest('Gamma');
            SingleTest('GenNor1');
            SingleTest('GenNor2');
            SingleTest('HyperbolicTan');
            SingleTest('JohnsonSB');
            SingleTest('JohnsonSU');
            SingleTest('Laplace');
            SingleTest('Logistic');
            SingleTest('Lognormal');
            SingleTest('LognormalMS');
            SingleTest('NakaRush');
            SingleTest('Normal');
            SingleTest('Pareto');
            SingleTest('Quantal');
            SingleTest('Quick');
            SingleTest('Rayleigh');
            SingleTest('Triangular');
            SingleTest('RNGamma');
            SingleTest('RNGammaMn');
            SingleTest('RNGammaMS');
            SingleTest('Rosin');
            SingleTest('t');
            SingleTest('TriangularCW');
            SingleTest('Uniform');
            SingleTest('UniformCW');
            SingleTest('UniGap');
            SingleTest('Wald2');
            SingleTest('Weibull');
            WriteMarker([MarkerRoot 'CntFas']);
        end
        
        if Continuous && Slow
            SingleTest('ChiSqNoncentral');
            SingleTest('ExpSumT');
            SingleTest('FNoncentral');
            SingleTest('r');
            SingleTest('rNoncentral');
            SingleTest('Recinormal');
            SingleTest('SkewNor');
            SingleTest('StudRng');
            SingleTest('tNoncentral');
            SingleTest('tPowerEst');
            SingleTest('TriangularG');
            SingleTest('TriangularGCWP');
            SingleTest('VonMises');
            SingleTest('Wald');
            WriteMarker([MarkerRoot 'CntSlo']);
        end
        
        if Discrete && Fast
            SingleTest('Binomial');
            SingleTest('BinomialMixed');
            SingleTest('Geometric');
            SingleTest('Poisson');
% NWJEFF SingleTest('RankDist');
            SingleTest('UniformInt');
            WriteMarker([MarkerRoot 'DisFas']);
        end
        
        if Discrete && Slow
            %    SingleTest('');
            WriteMarker([MarkerRoot 'DisSlo']);
        end
        
        if DerivedCont && Fast
            % This section takes ~14 min non-parallel
            SingleTest('AddTrans');
            SingleTest('ArcsinTrans');
            SingleTest('ExpTrans');
            SingleTest('InverseTrans');
            SingleTest('LinearTrans');
            SingleTest('LogisticTrans');
            SingleTest('LogitTrans');
            SingleTest('LogTrans');
            SingleTest('Mixture');
            SingleTest('MonotoneTrans');
            SingleTest('MultTrans');
            SingleTest('PhiTrans');
            SingleTest('PhiInvTrans');
            SingleTest('PowerTrans');
            SingleTest('Product');
            SingleTest('Ratio');
            SingleTest('SqrTrans');
            SingleTest('SqrtTrans');
            SingleTest('TruncatedP');
            SingleTest('TruncatedX');
            SingleTest('TruncatedXlow');
            SingleTest('TruncatedXhi');

            % Forward and reverse transformation combinations;  NWJEFF: Many more are possible.
            SingleTest('ExpLog');
            SingleTest('LogitLogistic');
            WriteMarker([MarkerRoot 'DerivContFas']);
        end
        
        if DerivedCont && Slow
            % This section takes 330 min non-parallel
            SingleTest('AttainP');
            SingleTest('Convolution');
            SingleTest('Difference');
            SingleTest('InfMix');
            SingleTest('MinBound');
            SingleTest('Order');
            SingleTest('OrderIID');
            WriteMarker([MarkerRoot 'DerivContSlo']);
        end
        
        if DerivedDisc && Fast
            % This section takes ~14 min non-parallel
            SingleTest('AddTransDisc');
            SingleTest('ArcsinTransDisc');
            SingleTest('ExpTransDisc');
            SingleTest('InverseTransDisc');
            SingleTest('LinearTransDisc');
            SingleTest('LogisticTransDisc');
            SingleTest('LogitTransDisc');
            SingleTest('LogTransDisc');
            SingleTest('MixtureDisc');
            SingleTest('MonotoneTransDisc');
            SingleTest('MultTransDisc');
            SingleTest('PhiInvTransDisc');
            SingleTest('PhiTransDisc');
            SingleTest('PowerTransDisc');
            SingleTest('ProductDisc');
            SingleTest('RatioDisc');
            SingleTest('SqrTransDisc');
            SingleTest('SqrtTransDisc');
            SingleTest('TruncatedPDisc');
            SingleTest('TruncatedXDisc');
            SingleTest('TruncatedXlowDisc');
            SingleTest('TruncatedXhiDisc');

            % Forward and reverse transformation combinations;  NWJEFF: Many more possible.
            SingleTest('ExpLogDisc');
            SingleTest('LogitLogisticDisc');
            WriteMarker([MarkerRoot 'DerivDiscFas']);
        end
        
        if DerivedDisc && Slow
%            SingleTest('AttainPDisc');
            SingleTest('ConvolutionDisc');
            SingleTest('DifferenceDisc');
%            SingleTest('InfMixDisc');
%            SingleTest('MinBoundDisc');
            SingleTest('OrderDisc');
            SingleTest('OrderIIDDisc');
            WriteMarker([MarkerRoot 'DerivDiscSlo']);
        end
        
        % Note: When running parallel, MATLAB reports xxx seconds testing time, but it appears to sum the
        % times for the different processors working in parallel.
        
        totalmins_elapsed = toc(time0) / 60
        
    end
    
    if (totalmins_elapsed > 5) && WantBeep
        DoneBeep;
    end
    
end  % function Run

function WriteMarker(fName)
f = fopen([fName '.chk'],'wt');
fclose(f);
end

