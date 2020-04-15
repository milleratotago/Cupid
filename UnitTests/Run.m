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
        
        if Continuous && Fast   % iCase 1
            rpt = cell(0,0);
            sectiontime = tic;
            rpt{end+1} = SingleTest('Beta');
            rpt{end+1} = SingleTest('Cauchy');
            rpt{end+1} = SingleTest('Chi');
            rpt{end+1} = SingleTest('ChiSq');
            rpt{end+1} = SingleTest('Cosine');
            rpt{end+1} = SingleTest('DblMon');
            rpt{end+1} = SingleTest('ExGauss');
            rpt{end+1} = SingleTest('ExGauMn');
            rpt{end+1} = SingleTest('ExGauRatio');
            rpt{end+1} = SingleTest('Exponential');
            rpt{end+1} = SingleTest('ExponenMn');
            rpt{end+1} = SingleTest('ExpSum');
            rpt{end+1} = SingleTest('ExtrVal1');
            rpt{end+1} = SingleTest('ExtrVal2');
            rpt{end+1} = SingleTest('ExtrValGen');
            rpt{end+1} = SingleTest('ExWald');
            rpt{end+1} = SingleTest('ExWaldMn');
            rpt{end+1} = SingleTest('ExWaldMSM');
            rpt{end+1} = SingleTest('F');
            rpt{end+1} = SingleTest('Gamma');
            rpt{end+1} = SingleTest('Geary');
            rpt{end+1} = SingleTest('GenNor1');
            rpt{end+1} = SingleTest('GenNor2');
            rpt{end+1} = SingleTest('HyperbolicTan');
            rpt{end+1} = SingleTest('JohnsonSB');
            rpt{end+1} = SingleTest('JohnsonSU');
            rpt{end+1} = SingleTest('KolmSmir');
            rpt{end+1} = SingleTest('Laplace');
            rpt{end+1} = SingleTest('Logistic');
            rpt{end+1} = SingleTest('Lognormal');
            rpt{end+1} = SingleTest('LognormalMCV');
            rpt{end+1} = SingleTest('LognormalMS');
            rpt{end+1} = SingleTest('NakaRush');
            rpt{end+1} = SingleTest('Normal');
            rpt{end+1} = SingleTest('Pareto');
            rpt{end+1} = SingleTest('ProdUni01');
            rpt{end+1} = SingleTest('Quantal');
            rpt{end+1} = SingleTest('Quick');
            rpt{end+1} = SingleTest('Rayleigh');
            rpt{end+1} = SingleTest('Triangular');
            rpt{end+1} = SingleTest('RNGamma');
            rpt{end+1} = SingleTest('RNGammaMn');
            rpt{end+1} = SingleTest('RNGammaMS');
            rpt{end+1} = SingleTest('Rosin');
            rpt{end+1} = SingleTest('t');
            rpt{end+1} = SingleTest('TriangularCW');
            rpt{end+1} = SingleTest('Uniform');
            rpt{end+1} = SingleTest('UniformCW');
            rpt{end+1} = SingleTest('UniGap');
            rpt{end+1} = SingleTest('Wald2');
            rpt{end+1} = SingleTest('Wald2MSD');
            rpt{end+1} = SingleTest('Weibull');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '1CntFas'],rpt,sectiontime);
        end
        
        if Discrete && Fast   % iCase 2
            rpt = cell(0,0);
            sectiontime = tic;
            rpt{end+1} = SingleTest('Binomial');
            rpt{end+1} = SingleTest('BinomialP');
            rpt{end+1} = SingleTest('BinomialZ');
            rpt{end+1} = SingleTest('BinomialMixed');
            rpt{end+1} = SingleTest('Geometric');
            rpt{end+1} = SingleTest('List');
            rpt{end+1} = SingleTest('NegativeBinomial');
            rpt{end+1} = SingleTest('Poisson');
            rpt{end+1} = SingleTest('RankDist');
            rpt{end+1} = SingleTest('UniformInt');
%            rpt{end+1} = SingleTest('YNdPrime');
%            rpt{end+1} = SingleTest('YNdPrimeSym');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '2DisFas'],rpt,sectiontime);
        end
        
        if DerivedCont && Fast   % iCase 3
            rpt = cell(0,0);
            sectiontime = tic;
            rpt{end+1} = SingleTest('AddTrans');
            rpt{end+1} = SingleTest('ArcsinTrans');
            rpt{end+1} = SingleTest('ExpTrans');
            rpt{end+1} = SingleTest('InverseTrans');
            rpt{end+1} = SingleTest('LinearTrans');
            rpt{end+1} = SingleTest('LikeRatioC');
            rpt{end+1} = SingleTest('LnLikeRatioC');
            rpt{end+1} = SingleTest('LogisticTrans');
            rpt{end+1} = SingleTest('LogitTrans');
            rpt{end+1} = SingleTest('LogTrans');
            rpt{end+1} = SingleTest('Mixture');
            rpt{end+1} = SingleTest('MonotoneTrans');
            rpt{end+1} = SingleTest('MultTrans');
            rpt{end+1} = SingleTest('PhiTrans');
            rpt{end+1} = SingleTest('PhiInvTrans');
            rpt{end+1} = SingleTest('PowerTrans');
            rpt{end+1} = SingleTest('Product');
            rpt{end+1} = SingleTest('ProdUni0p');
            rpt{end+1} = SingleTest('Ratio');
            rpt{end+1} = SingleTest('SqrTrans');
            rpt{end+1} = SingleTest('SqrtTrans');
            rpt{end+1} = SingleTest('TruncatedP');
            rpt{end+1} = SingleTest('TruncatedX');
            rpt{end+1} = SingleTest('TruncatedXlow');
            rpt{end+1} = SingleTest('TruncatedXhi');

            % Forward and reverse transformation combinations;  NWJEFF: Many more are possible.
            rpt{end+1} = SingleTest('ExpLog');
            rpt{end+1} = SingleTest('LogitLogistic');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '3DerivContFas'],rpt,sectiontime);
        end
        
        if DerivedDisc && Fast   % iCase 4
            rpt = cell(0,0);
            sectiontime = tic;
            rpt{end+1} = SingleTest('AddTransDisc');
            rpt{end+1} = SingleTest('ArcsinTransDisc');
            rpt{end+1} = SingleTest('ExpTransDisc');
            rpt{end+1} = SingleTest('InverseTransDisc');
            rpt{end+1} = SingleTest('LikeRatioD');
            rpt{end+1} = SingleTest('LnLikeRatioD');
            rpt{end+1} = SingleTest('LinearTransDisc');
            rpt{end+1} = SingleTest('LogisticTransDisc');
            rpt{end+1} = SingleTest('LogitTransDisc');
            rpt{end+1} = SingleTest('LogTransDisc');
            rpt{end+1} = SingleTest('MixtureDisc');
            rpt{end+1} = SingleTest('MonotoneTransDisc');
            rpt{end+1} = SingleTest('MultTransDisc');
            rpt{end+1} = SingleTest('PhiInvTransDisc');
            rpt{end+1} = SingleTest('PhiTransDisc');
            rpt{end+1} = SingleTest('PowerTransDisc');
            rpt{end+1} = SingleTest('ProductDisc');
            rpt{end+1} = SingleTest('RatioDisc');
            rpt{end+1} = SingleTest('SqrTransDisc');
            rpt{end+1} = SingleTest('SqrtTransDisc');
            rpt{end+1} = SingleTest('TruncatedPDisc');
            rpt{end+1} = SingleTest('TruncatedXDisc');
            rpt{end+1} = SingleTest('TruncatedXlowDisc');
            rpt{end+1} = SingleTest('TruncatedXhiDisc');

            % Forward and reverse transformation combinations;  NWJEFF: Many more possible.
            rpt{end+1} = SingleTest('ExpLogDisc');
            rpt{end+1} = SingleTest('LogitLogisticDisc');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '4DerivDiscFas'],rpt,sectiontime);
        end
        
        if Continuous && Slow   % iCase 5
            rpt = cell(0,0);
            sectiontime = tic;
            rpt{end+1} = SingleTest('ChiSqNoncentral');
            rpt{end+1} = SingleTest('dMATLABc');
            rpt{end+1} = SingleTest('ExpSumT');
            rpt{end+1} = SingleTest('FNoncentral');
            rpt{end+1} = SingleTest('ghHoag');
            rpt{end+1} = SingleTest('r');
            rpt{end+1} = SingleTest('rNoncentral');
            rpt{end+1} = SingleTest('Recinormal');
            rpt{end+1} = SingleTest('RecinormalInv');
            rpt{end+1} = SingleTest('RecinormalMS');
            rpt{end+1} = SingleTest('SkewNor');
            rpt{end+1} = SingleTest('StochCasc2T');
            rpt{end+1} = SingleTest('StudRng');
            rpt{end+1} = SingleTest('tNoncentral');
            rpt{end+1} = SingleTest('tPowerEst');
            rpt{end+1} = SingleTest('TriangularG');
            rpt{end+1} = SingleTest('TriangularGCWP');
            rpt{end+1} = SingleTest('VonMises');
            rpt{end+1} = SingleTest('Wald');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '5CntSlo'],rpt,sectiontime);
        end
        

        if Discrete && Slow   % iCase 6
            rpt = cell(0,0);
            sectiontime = tic;
            %    rpt{end+1} = SingleTest('');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '6DisSlo'],rpt,sectiontime);
        end
        

        if DerivedCont && Slow   % iCase 7
            rpt = cell(0,0);
            sectiontime = tic;
            rpt{end+1} = SingleTest('AttainP');
            rpt{end+1} = SingleTest('ConditXLTY');
            rpt{end+1} = SingleTest('ConditXGTY');
            rpt{end+1} = SingleTest('Convolution');
            rpt{end+1} = SingleTest('Difference');
            rpt{end+1} = SingleTest('InfMix');
            rpt{end+1} = SingleTest('MinBound');
            rpt{end+1} = SingleTest('Order');
            rpt{end+1} = SingleTest('OrderIID');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '7DerivContSlo'],rpt,sectiontime);
        end
        

        if DerivedDisc && Slow   % iCase 8
            rpt = cell(0,0);
            sectiontime = tic;
%            rpt{end+1} = SingleTest('AttainPDisc');
            rpt{end+1} = SingleTest('ConvolutionDisc');
            rpt{end+1} = SingleTest('DifferenceDisc');
%            rpt{end+1} = SingleTest('InfMixDisc');
%            rpt{end+1} = SingleTest('MinBoundDisc');
            rpt{end+1} = SingleTest('OrderDisc');
            rpt{end+1} = SingleTest('OrderIIDDisc');
            sectiontime = toc(sectiontime) / 60;
            WriteMarker([MarkerRoot '8DerivDiscSlo'],rpt,sectiontime);
        end
        
        % Note: When running parallel, MATLAB reports xxx seconds testing time, but it appears to sum the
        % times for the different processors working in parallel.
        
        totalmins_elapsed = toc(time0) / 60
        
    end
    
    if (totalmins_elapsed > 5) && WantBeep
        DoneBeep;
    end
    
end  % function Run

function WriteMarker(fName,rpt,min)
f = fopen([fName '.chk'],'wt');
fprintf(f,'%s\n',rpt{:});
fprintf(f,'%10.3f total minutes for this section.\n',min);
fclose(f);
end

