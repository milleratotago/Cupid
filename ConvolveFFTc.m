classdef ConvolveFFTc < dContinuous
    % ConvolveFFTc(varargin) creates a random variable that is an FFT-based approximation
    %  of the sum of k independent, continuous basis random variables.
    % NEWJEFF: Right now this uses spline interpolation but linear interp1 might be better.
    
    % Optional parameters:
    %   NxPoints: Number of points used to approximate the full range of X
    %   delx:     The resolution used to approximate the full range of X
    %     Notes: If both of these 2 optional parameters are specified, delx is used.
    %            If neither of these 2 optional parameters are specified, defaultNxPoints is used.
    % SplinePDF is used.
    
    properties (Constant)
        defaultNxPoints = 2^12;  % 4096 is a reasonable speed/accuracy compromise in cases examined so far.
    end
    
    properties(SetAccess = public)
        NxPoints
    end

    properties(SetAccess = private)
        delx
        BasisRV
        NDists
        BasisMeans
        BasisVariances
        KnownMeans
        KnownVariances
        % The following are used for a discrete approximation
        % instead of Splines
        DiscreteXs
        DiscretePDFs
        DiscreteCDFs
    end
    
    methods
        
        function obj=ConvolveFFTc(varargin)
            obj=obj@dContinuous('ConvolveFFTc');
            varargin = obj.ExtractNxPoints(varargin{:});
            obj.BasisRV = varargin;
            obj.NDists = numel(obj.BasisRV);
            obj.NDistParms = 0;
            obj.DefaultParmCodes = '';
            for iDist=1:obj.NDists
                obj.DefaultParmCodes = [obj.DefaultParmCodes obj.BasisRV{iDist}.DefaultParmCodes];
                obj.NDistParms = obj.NDistParms + obj.BasisRV{iDist}.NDistParms;
            end
            obj.BasisMeans = zeros(obj.NDists,1);
            obj.BasisVariances = zeros(obj.NDists,1);
            obj.KnownMeans = false;
            obj.KnownVariances = false;
            obj.ReInit;
        end
        
        function varargin=ExtractNxPoints(obj,varargin)
            [obj.NxPoints, varargin] = ExtractNameVali('NxPoints',-1,varargin);
            [obj.delx, varargin] = ExtractNameVali('delx',-1,varargin);
            if (obj.NxPoints<0) &&  (obj.NxPoints<0)
                obj.NxPoints = ConvolveFFTc.defaultNxPoints;
            end
        end
        
        function BuildMyName(obj)
            sTemp = obj.BasisRV{1}.StringName;
            for iDist=2:obj.NDists
                sTemp = [sTemp ',' obj.BasisRV{iDist}.StringName]; %#ok<AGROW>
            end
            obj.StringName = [obj.FamilyName '(' sTemp ')'];
        end
        
        function []=ClearBeforeResetParmsC(obj)
            ClearBeforeResetParmsC@dContinuous(obj);
            obj.KnownMeans = false;
            obj.KnownVariances = false;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            ParmsUsed = 0;
            for iDist = 1:obj.NDists
                obj.BasisRV{iDist}.ResetParms(newparmvalues(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms));
                if obj.NameBuilding
                    obj.BasisRV{iDist}.BuildMyName;
                end
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            ParmsUsed = 0;
            for iDist = 1:obj.NDists
                obj.BasisRV{iDist}.PerturbParms(ParmCodes(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms));
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
            obj.ReInit;
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [];
            ParmsUsed = 0;
            for iDist=1:obj.NDists
                Reals = [Reals obj.BasisRV{iDist}.ParmsToReals(Parms(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms))]; %#ok<AGROW>
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [];
            ParmsUsed = 0;
            for iDist=1:obj.NDists
                Parms = [Parms obj.BasisRV{iDist}.RealsToParms(Reals(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms))]; %#ok<AGROW>
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [];
            for iDist=1:obj.NDists
                parmvals = [parmvals obj.BasisRV{iDist}.ParmValues]; %#ok<AGROW>
            end
        end
        
        function []=ReInit(obj)
            % Using existing bases, find bounds, set up Xs, run FFT
            % 2022-07-06 bombs in this routine because of empty pos_pdfs associated with SD=0 for one basis.
            % 2022-12-06 bombs for exGamma(2.568147394198137e-07,63.335183209196250,625.1190965420082) because of empty pos_pdfs
            
            % collect lower/upper bound info about bases:
            lowerb = zeros(obj.NDists,1);
            upperb = zeros(obj.NDists,1);
            for i=1:obj.NDists
                lowerb(i) = obj.BasisRV{i}.LowerBound;
                upperb(i) = obj.BasisRV{i}.UpperBound;
            end
            % Compute corresponding properties of convolution
            lowerconv = sum(lowerb);
            upperconv = sum(upperb);
            
            % Compute full range of X's spanning all bases & convolution
            % The X's will reflect shifted versions of all RVs
            % (i.e., shifted to have lower bounds of 0.
            x_min = 0;
            x_max = max([sum(upperb-lowerb); upperconv-lowerconv]);
            
            % Compute the actual X's
            if obj.delx > 0
                useNxPoints = ceil( (x_max - x_min) / obj.delx );
            else
                useNxPoints = obj.NxPoints;
            end
            X = linspace(x_min,x_max,useNxPoints);
            
            % 2022-07-07 experimental NEWJEFF
            % if there are only 2 distributions and one of them accounts for >99% of the variance,
            % just use the PDFs for that one (adjusting for the constant shift of the other).
            if obj.NDists == 2
                Var1 = obj.BasisRV{1}.Variance;
                Var2 = obj.BasisRV{2}.Variance;
                TtlVar = Var1 + Var2;
                MostlyVar1 = Var1 / TtlVar > 0.99;
                MostlyVar2 = Var2 / TtlVar > 0.99;
            else
                MostlyVar1 = false;
                MostlyVar2 = false;
            end
            if MostlyVar1
                ConstMu = obj.BasisRV{2}.Mean;
                obj.DiscretePDFs = obj.BasisRV{1}.PDF(X-ConstMu);
            end
            if MostlyVar2
                ConstMu = obj.BasisRV{1}.Mean;
                obj.DiscretePDFs = obj.BasisRV{2}.PDF(X-ConstMu);
            end
            if MostlyVar1 || MostlyVar2
                obj.DiscreteXs = X;
                obj.DiscreteCDFs = cumsum(obj.DiscretePDFs);
                obj.DiscreteCDFs = obj.DiscreteCDFs / obj.DiscreteCDFs(end);
                xStep = obj.DiscreteXs(2) - obj.DiscreteXs(1)           ;
                obj.DiscretePDFs = obj.DiscretePDFs / xStep;
                obj.LowerBound = X(1);
                obj.UpperBound = X(end);
                obj.Initialized = true;
                if obj.NameBuilding
                    obj.BuildMyName;
                end
                return
            end
            
            % Load Y with PDFs of bases:
            % I don't understand the various discussions about using fftshift & ifftshift
            % but it seems correct without them.
            Y = zeros(obj.NDists,useNxPoints);
            for i=1:obj.NDists
                % Add lowerb for pdf computation to compensate for
                % the use of shifted X values
                % Y(i,:) = ifftshift(obj.BasisRV{i}.PDF(X+lowerb(i)));
                Y(i,:) = obj.BasisRV{i}.PDF(X+lowerb(i));
            end
            Y(isnan(Y)) = 0;  % ifft produces all nan's if any Y is nan
            
%             % NEWJEFF: Experimental
%             % Get rid of the Y's for any RV with too-small variance
%             % because these cause problems for fft
%             totalVar = 0;
%             for iDist=1:obj.NDists
%                 totalVar = totalVar + obj.BasisRV{iDist}.Variance;
%             end
%             tooSmallVar = false(obj.NDists,1);
%             for iDist=1:obj.NDists
%                 tooSmallVar(iDist) = obj.BasisRV{iDist}.SD < 0.001*sqrt(totalVar);
%             end
%             Y(tooSmallVar,:) = [];
            
            % Compute PDF of convolution via FFT
            % I don't really understand this
            fftpdf = real(ifft(prod(fft(Y,[],2))));  % Problem: Could be all zeros and nans
            % fftpdf = fftshift(fftpdf);
            % NEWJEFF: Following line seems redundant with later normalization
            fftpdf = fftpdf / trapz(X,fftpdf);  % normalize pdf values for integral=1
            
            % At this point the fftpdf values are correct, but the X vector still
            % reflects the shifts, so remove those.
            needed_shift = lowerconv;
            X = X + needed_shift;
            
            % Trim X and fftpdf vectors to get rid of leading & trailing zeros, if any.
            % (zero's at each end are added below)
            pos_pdfs = find(fftpdf>0);
            % Bombs here if fftpdf has no entries greater than 0 (presumably all NaN)
            first_pos_pdf = pos_pdfs(1);
            last_pos_pdf = pos_pdfs(end);
            X = X(first_pos_pdf:last_pos_pdf);
            fftpdf = fftpdf(first_pos_pdf:last_pos_pdf);
            
%             % These are the final X and fftpdf values used for splinepdf
%             obj.SplinePDFsXs = X;
%             obj.SplinePDFs = fftpdf;
%             obj.PDFSplineInfo = spline(obj.SplinePDFsXs,[0 obj.SplinePDFs 0]);  % Include end 0's to guarantee flatness at edges.
%             obj.UseSplinePDF = true;
%             obj.HaveSplinePDFs = true;
            
            % These are the final X, pdf, and cdf values used for the discrete approximation:
            fftpdf(fftpdf<0) = 0;  % Must be numerical errors so get rid of them
            fftpdf = fftpdf / sum(fftpdf); %  / trapz(X,fftpdf);  % normalize to area 1
            obj.DiscreteXs = X;
            obj.DiscretePDFs = fftpdf;
            obj.DiscreteCDFs = cumsum(obj.DiscretePDFs);
            obj.DiscreteCDFs = obj.DiscreteCDFs / obj.DiscreteCDFs(end);
            xStep = obj.DiscreteXs(2) - obj.DiscreteXs(1)           ;
            obj.DiscretePDFs = obj.DiscretePDFs / xStep;
            
            obj.LowerBound = X(1);
            obj.UpperBound = X(end);
            obj.Initialized = true;
            if obj.NameBuilding
                obj.BuildMyName;
            end
        end
        
        function thispdf = PDF(obj,X)
            thispdf=zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thispdf(InBounds) = interp1(obj.DiscreteXs,obj.DiscretePDFs,X(InBounds));
        end
        
        function thiscdf = CDF(obj,X)
            thiscdf = zeros(size(X));
            thiscdf(X>=obj.UpperBound) = 1;
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(InBounds) = interp1(obj.DiscreteXs,obj.DiscreteCDFs,X(InBounds));
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.KnownMeans
                for iDist=1:obj.NDists
                    obj.BasisMeans(iDist) = obj.BasisRV{iDist}.Mean;
                end
                obj.KnownMeans = true;
            end
            thisval = sum(obj.BasisMeans);
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.KnownVariances
                for iDist=1:obj.NDists
                    obj.BasisVariances(iDist) = obj.BasisRV{iDist}.Variance;
                end
                obj.KnownVariances = true;
            end
            thisval = sum(obj.BasisVariances);
        end
    
        function thisVal = Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            basisVals = zeros(varargin{:});  % make an array of the required size for each basis
            nPerBasis = numel(basisVals);
            AllBases = zeros(obj.NDists,nPerBasis);
            for iDist=1:obj.NDists
                AllBases(iDist,:) = obj.BasisRV{iDist}.Random(1,nPerBasis);
            end
            sums = sum(AllBases,1);
            thisVal = reshape(sums,size(basisVals));
        end
        
end  % methods

end  % class ConvolveFFTc

