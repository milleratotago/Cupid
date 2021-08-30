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
        defaultNxPoints = 2^16;
    end
    
    properties(SetAccess = private)
        NxPoints
        delx
        BasisRV
        NDists
        BasisMeans
        BasisVariances
        KnownMeans
        KnownVariances
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
        
        function []=ClearBeforeResetParms(obj)
            ClearBeforeResetParms@dContinuous;
            obj.KnownMeans = false;
            obj.KnownVariances = false;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
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
            
            % collect lower/upper bound & mean info about bases:
            lowerb = zeros(obj.NDists,1);
            upperb = zeros(obj.NDists,1);
            for i=1:obj.NDists
                lowerb(i) = obj.BasisRV{i}.LowerBound;
                upperb(i) = obj.BasisRV{i}.UpperBound;
            end
            % Compute corresponding properties of convolution
            lowerconv = sum(lowerb);
            upperconv = sum(upperb);
            %             meanconv = sum(obj.BasisMeans);
            
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
            
            % Compute PDF of convolution via FFT
            % I don't really understand this
            fftpdf = real(ifft(prod(fft(Y,[],2))));
            % fftpdf = fftshift(fftpdf);
            fftpdf = fftpdf / trapz(X,fftpdf);  % normalize pdf values for integral=1
            
            % At this point the fftpdf values are correct, but the X vector still
            % reflects the shifts, so remove those.
            needed_shift = sum(lowerb);
            X = X + needed_shift;
            
            % Trim X and fftpdf vectors to get rid of leading & trailing zeros, if any.
            % (zero's at each end are added below)
            pos_pdfs = find(fftpdf>0);
            first_pos_pdf = pos_pdfs(1);
            last_pos_pdf = pos_pdfs(end);
            X = X(first_pos_pdf:last_pos_pdf);
            fftpdf = fftpdf(first_pos_pdf:last_pos_pdf);
            
            % These are the final X and fftpdf values used for splinepdf
            obj.SplinePDFsXs = X;
            obj.SplinePDFs = fftpdf;
            obj.PDFSplineInfo = spline(obj.SplinePDFsXs,[0 obj.SplinePDFs 0]);  % Include end 0's to guarantee flatness at edges.
            obj.UseSplinePDF = true;
            obj.HaveSplinePDFs = true;
            
            obj.LowerBound = X(1);
            obj.UpperBound = X(end);
            obj.Initialized = true;
            if obj.NameBuilding
                obj.BuildMyName;
            end
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

