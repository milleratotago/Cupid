classdef ConvolveFFTcIID < dContinuous
    % ConvolveFFTcIID(NDists,Basis,varargin)
    % Convolution of NDists i.i.d. random variables with distribution BasisRV{1}
    % NEWJEFF: There is a lot of code overlap with ConvolveFFTc, but it will be tough
    % to do this right (even just initial constructors!).
    % NEWJEFF: It would be nice to add a discrete approx. for convolveFFTcIID like convolveFFTc

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
        
        function obj=ConvolveFFTcIID(NDists,Basis,varargin)
            obj=obj@dContinuous('ConvolveFFTcIID');
            obj.ExtractNxPoints(varargin{:});
            obj.BasisRV = Basis;
            obj.NDists = NDists;
            obj.NDistParms = obj.BasisRV.NDistParms;
            obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
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
            obj.StringName = [obj.FamilyName '(' num2str(obj.NDists) ',' obj.BasisRV.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.BasisRV.ResetParms(newparmvalues);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            obj.ReInit;
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = obj.BasisRV.ParmsToReals(Parms);
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = obj.BasisRV.RealsToParms(Reals);
        end
        
        function parmvals = ParmValues(obj)
            parmvals = obj.BasisRV.ParmValues;
        end
        
        function []=ReInit(obj)
            % Using existing bases, find bounds, set up Xs, run FFT
            
            % get lower/upper bound info about basis:
            lowerb = obj.BasisRV.LowerBound;
            upperb = obj.BasisRV.UpperBound;
            % Compute corresponding properties of convolution
            lowerconv = lowerb * obj.NDists;
            upperconv = upperb * obj.NDists;
            
            % Compute full range of X's spanning all bases & convolution
            % The X's will reflect shifted versions of all RVs
            % (i.e., shifted to have lower bounds of 0.
            x_min = 0;
            x_max = max([upperb-lowerb; upperconv-lowerconv]);
            
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
                % Y(i,:) = ifftshift(obj.BasisRV.PDF(X+lowerb));
                if i==1
                    Y(i,:) = obj.BasisRV.PDF(X+lowerb);
                 else
                    Y(i,:) = Y(i-1,:);
                 end
            end
            
            % Compute PDF of convolution via FFT
            % I don't really understand this
            fftpdf = real(ifft(prod(fft(Y,[],2))));
            % fftpdf = fftshift(fftpdf);
            fftpdf = fftpdf / trapz(X,fftpdf);  % normalize pdf values for integral=1
            
            % At this point the fftpdf values are correct, but the X vector still
            % reflects the shifts, so remove those.
            needed_shift = lowerconv;
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
                obj.BasisMeans = obj.BasisRV.Mean;
                obj.KnownMeans = true;
            end
            thisval = obj.NDists * obj.BasisMeans;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if ~obj.KnownVariances
                obj.BasisVariances = obj.BasisRV.Variance;
                obj.KnownVariances = true;
            end
            thisval = obj.NDists * obj.BasisVariances;
        end
    
        function thisVal = Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            basisVals = zeros(varargin{:});  % make an array of the required size for each basis
            nPerBasis = numel(basisVals);
            AllBases = zeros(obj.NDists,nPerBasis);
            for iDist=1:obj.NDists
                AllBases(iDist,:) = obj.BasisRV.Random(1,nPerBasis);
            end
            sums = sum(AllBases,1);
            thisVal = reshape(sums,size(basisVals));
        end
        
    end % methods

end

