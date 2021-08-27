classdef utdMATLABc < utContinuous
    % Unit testing for the dMATLABc class, which allows continuous MATLAB probability distribution objects
    % to be used as Cupid distributions.
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23};
        % parmCase = {11};  % Check a single distribution
        % Special handling for some distributions:
        %  8: GeneralizedPareto: Do not adjust constant to avoid impossible data values.
        %  9: HalfNormal: Do not adjust mu to avoid impossible data values.
        % 18: Rician; Modified ErrFound to ignore many complaints that NCX2INV did not converge
        % 19: Stable: Moments only exist with alpha=2, & estimation problems when beta is free to vary.
        % 20: tLocationScale: nu held constant
        % 21: Triangular: Skip MLEst to avoid checking impossible data values
        % 22: Uniform: Skip MLEst to avoid checking impossible data values.
        % Also: Triangular and Uniform parameter constraints are ignored.
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utdMATLABc(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the dMATLABc distribution.

            testCase.ThisCase  = parmCase;
            % Initialize distributions/parameters & adjust tolerances as appropriate:
            switch parmCase
                case 1
                    pd = makedist('Beta','a',2,'b',3);
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.010);
                case 2
                    pd = makedist('BirnbaumSaunders');
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                case 3
                    pd = makedist('Burr','a',1,'c',9,'k',2);
                    testCase.Dist = dMATLABc(pd,'rrr',[0 0 0],[+inf +inf +inf]);
                    SetTolerances(testCase,0.010);
                    testCase.ParmEstAbsTol = 0.025;  % Parameter estimation is not so accurate.
                    testCase.ParmEstRelTol = 0.025;
                case 4
                    pd = makedist('Exponential','mu',100);
                    testCase.Dist = dMATLABc(pd,'r',0,+inf);
                    SetTolerances(testCase,0.005);
                    testCase.ParmEstAbsTol = 0.01;  % Parameter estimation is not so accurate.
                    testCase.ParmEstRelTol = 0.01;
                case 5
                    pd = makedist('ExtremeValue');
                    testCase.Dist = dMATLABc(pd,'rr',[-inf 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                case 6
                    pd = makedist('Gamma');
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.010);
                case 7
                    pd = makedist('GeneralizedExtremeValue');
                    testCase.Dist = dMATLABc(pd,'rrr',[-inf 0 -inf],[+inf +inf +inf]);
                    SetTolerances(testCase,0.010);
                case 8
                    pd = makedist('GeneralizedPareto','k',0.2,'sigma',1,'theta',0);
                    testCase.Dist = dMATLABc(pd,'rrf',[-inf 0 -inf],[+inf +inf +inf]);
                    % Do not adjust offset parameter theta during estimation,
                    % because when it is adjusted some data values are impossible.
                    SetTolerances(testCase,0.007);
                case 9
                    pd = makedist('HalfNormal');
                    testCase.Dist = dMATLABc(pd,'fr',[-inf 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                    % testCase.SkipEstAll = true;  % When bounds are adjusted some data values are impossible.
                case 10
                    pd = makedist('InverseGaussian');
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.008);
                case 11
                    pd = makedist('Logistic','mu',10,'sigma',1);
                    testCase.Dist = dMATLABc(pd,'rr',[-inf 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                case 12
                    pd = makedist('Loglogistic','mu',2,'sigma',.2);
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                case 13
                    pd = makedist('Lognormal');
                    testCase.Dist = dMATLABc(pd,'rr',[-inf 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                    testCase.ParmEstAbsTol = 0.01;  % Parameter estimation is not so accurate.
                    testCase.ParmEstRelTol = 0.01;
                case 14
                    pd = makedist('Nakagami');
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.007);
                    testCase.ParmEstAbsTol = 0.01;  % Parameter estimation is not so accurate.
                    testCase.ParmEstRelTol = 0.01;
                case 15
                    pd = makedist('Normal','mu',0,'sigma',1);
                    testCase.Dist = dMATLABc(pd,'rr',[-inf 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                case 16
                    x = 0:10;
                    Fx = x.^2/100;
                    pd = makedist('PiecewiseLinear','x',x,'Fx',Fx);
                    testCase.Dist = dMATLABc(pd); % ,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                    testCase.SkipEstAll = true;  % There are no parameters
                case 17
                    pd = makedist('Rayleigh');
                    testCase.Dist = dMATLABc(pd,'r',0,+inf);
                    SetTolerances(testCase,0.005);
                case 18
                    pd = makedist('Rician','s',8,'sigma',1);
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[inf inf]);
                    SetTolerances(testCase,0.010);
                case 19
                    pd = makedist('Stable','alpha',2,'beta',0.25,'gam',1,'delta',0);    % Moments only exist when alpha=2
                    testCase.Dist = dMATLABc(pd,'ffrr',[0 -1 0 -inf],[2 1 +inf +inf]);  % Fix some parameters to speed searches & avoid estimation problems.
                    SetTolerances(testCase,0.010);
                case 20
                    pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',25);
                    testCase.Dist = dMATLABc(pd,'rrf',[-inf 0 0],[+inf +inf +inf]);
                    SetTolerances(testCase,0.010);
                case 21
                    pd = makedist('Triangular'); % NEWJEFF: Ignoring parameter constraints: lower limit, peak, upper limit
                    testCase.Dist = dMATLABc(pd,'rrr',[-inf -inf -inf],[+inf +inf +inf]);
                    SetTolerances(testCase,0.002);
                    testCase.SkipEstML = true;  % When bounds are adjusted some data values are impossible.
                case 22
                    pd = makedist('Uniform');  % NEWJEFF: Ignoring parameter constraints
                    testCase.Dist = dMATLABc(pd,'rr',[-inf -inf],[+inf +inf]);
                    SetTolerances(testCase,0.002);
                    testCase.SkipEstML = true;  % When bounds are adjusted some data values are impossible.
                case 23
                    pd = makedist('Weibull');
                    testCase.Dist = dMATLABc(pd,'rr',[0 0],[+inf +inf]);
                    SetTolerances(testCase,0.005);
                    testCase.ParmEstAbsTol = 0.01;  % Parameter estimation is not so accurate.
                    testCase.ParmEstRelTol = 0.01;
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            % SetupXs(testCase,40,200);
            testCase.xvalues = testCase.Dist.XsToPlot;
            testCase.xMLE = testCase.Dist.Random(100000,1);

            % testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            % testCase.Dist.SearchOptions.MaxIter = 20000;
            % testCase.RawMomentAbsTol(4)=.005; % Skewness computations slightly problematic
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
%         function ComparePDFs(testCase)
%             % Check matches to known PDFs
%             if testCase.ThisCase == 1
%                 tempU = Uniform(testCase.Dist.LowerBound,testCase.Dist.UpperBound);
%                 tempPDF = tempU.PDF(testCase.xvalues);
%                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
%             elseif testCase.ThisCase == 2
%                 tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
%                 tempPDF = tempU.PDF(testCase.xvalues);
%                 testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
%             elseif testCase.ThisCase == 3
%             end
%         end
        
    end
    
end  % utdMATLABc


