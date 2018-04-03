classdef ArcsinTrans < dMonoTrans0
    % ArcsinTrans(BasisRV): Arcsin(sqrt( Standard_Normal.CDF(BasisRV/2) ))
    % This transformation is used with frequency data (e.g. binomial proportions)
    % in order to obtain more normally distributed values.  For example,
    % suppose some p's have the skewed distribution b1=Beta(21,2), with RawSkewness = -0.06
    % The distribution ArcsinTrans(b1) is less skewed, RawSkewness = -0.012
    % (though RelSkewness is actually slightly larger: -1.149 vs -1.166).
    
    properties(SetAccess = protected)
        Standard_Normal
    end
    
    methods
        
        function obj=ArcsinTrans(BasisDist)
            obj=obj@dMonoTrans0('ArcsinTrans',BasisDist);
            obj.Standard_Normal = Normal(0,1);
            obj.TransReverses = false;
            obj.ReInit;
        end
        
        function Trans = PreTransToTrans(obj,PreTrans)
            % The "arcsin transformation" of a random variable.
            NCDF = obj.Standard_Normal.CDF(PreTrans/2);
            NCDF = sqrt(NCDF);
            Trans = asin(NCDF);
        end
        
        function PreTrans = TransToPreTrans(obj,Trans)
            % The "inverse arcsin transformation" of a random variable.
            Temp = sin(Trans);
            Temp = Temp.^2;
            PreTrans = 2*obj.Standard_Normal.InverseCDF(Temp);
        end
        
    end  % methods
    
end  % class ArcsinTrans


