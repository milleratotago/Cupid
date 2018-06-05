classdef RankDist < dDiscrete
    % RankDist(SampleSize,BasisX,BasisY): The probability distribution of the rank of X within a sample consisting
    % of one value of X and (SampleSize-1) values of Y.  RANK 1 MEANS LARGEST!
    
    % Notes:
    % NewJeff: Bad numerical errors lead to negative PDFs with RankDist(30,ChiSq(80),ChiSq(40))
    % NewJeff: This version only works with continuous Basis distributions.  With discrete ones, ties would be a problem.
    % o By default SampleSize is NOT adjusted when parameter fitting.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        SampleSize, BasisX, BasisY
    end
    
    methods
        
        function obj=RankDist(varargin)   % Constructor
            obj=obj@dDiscrete('RankDist');
            switch nargin
                case 0
                case 3
                    obj.BasisX = varargin{2};
                    obj.BasisY = varargin{3};
                    ResetParms(obj,[varargin{1} obj.BasisX.ParmValues obj.BasisY.ParmValues]);
                otherwise
                    ME = MException('RankDist:Constructor', ...
                        'RankDist constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % SampleSize is fixed.
            obj.BasisX.PerturbParms(ParmCodes(2:obj.BasisX.NDistParms+1));
            obj.BasisY.PerturbParms(ParmCodes(obj.BasisX.NDistParms+2:end));
            obj.ResetParms([obj.SampleSize obj.BasisX.ParmValues obj.BasisY.ParmValues]);
        end
        
        function BuildMyName(obj)
            obj.StringName = ['RankDist(' num2str(obj.SampleSize) ',' obj.BasisX.StringName  ',' obj.BasisY.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Initialized = false;
            obj.SampleSize = round(newparmvalues(1));  % Make sure it is an integer because fminsearch passes a real
            obj.BasisX.ResetParms(newparmvalues(2:obj.BasisX.NDistParms+1));
            obj.BasisY.ResetParms(newparmvalues(obj.BasisX.NDistParms+2:end));
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            obj.NValues = obj.SampleSize;
            obj.NDistParms = 1 + obj.BasisX.NDistParms + obj.BasisY.NDistParms;
            obj.DefaultParmCodes = ['f' obj.BasisX.DefaultParmCodes obj.BasisY.DefaultParmCodes];
            obj.LowerBound = 1;
            obj.UpperBound = obj.SampleSize;
            obj.MakeTables;
            obj.TrimTables(obj.PDFNearlyZero,1);
            obj.SetBinEdges;
            obj.Initialized = true;
            if obj.NameBuilding
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(2,Parms(1)) obj.BasisX.ParmsToReals(Parms(2:obj.BasisX.NDistParms+1)) obj.BasisY.ParmsToReals(Parms(obj.BasisX.NDistParms+2:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(2,Reals(1)) obj.BasisX.RealsToParms(Reals(2:obj.BasisX.NDistParms+1)) obj.BasisY.RealsToParms(Reals(obj.BasisX.NDistParms+2:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.SampleSize obj.BasisX.ParmValues obj.BasisY.ParmValues];
        end
        
        function []=MakeTables(obj)
            m = obj.SampleSize;  % Make sure it is an integer because fminsearch passes a real
            
            % Compute prob X has rank one in subsamples with 0..m-1 Y's.
            PrXisMaxInSubset = ones(m,1);
            PrXisMaxInSubset(2) = PrXGTY(obj.BasisX,obj.BasisY);
            for k=3:m
                MaxOfSubset = OrderIID(k-1,k-1,obj.BasisY);
                PrXisMaxInSubset(k) = PrXGTY(obj.BasisX,MaxOfSubset);
            end
            
            obj.DiscreteX = 1:m;
            obj.DiscretePDF = zeros(1,m);
            sump = 0;
            for k = 1:m
                thisx = k;
                obj.DiscreteX(k) = thisx;
                thispdf = 0;
                for i=0:k-1
                    thispdf = thispdf + nchoosek(k-1,i)*(-1)^i*PrXisMaxInSubset(m-k+i+1);
                end
                obj.DiscretePDF(k) = nchoosek(m-1,k-1)*thispdf;
            end
            obj.DiscreteCDF = cumsum(obj.DiscretePDF);
            obj.DiscreteCDF(end) = 1;
            obj.StoredTablesInitialized = true;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval=zeros(varargin{:});
            for i=1:numel(thisval)
                xsample = Random(obj.BasisX);
                ysample = Random(obj.BasisY,obj.SampleSize-1,1);
                thisval(i) = obj.SampleSize-sum(xsample>ysample);
            end
        end
        
    end  % methods
    
end  % class RankDist

