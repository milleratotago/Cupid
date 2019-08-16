classdef KolmSmir < dContinuous
    % Kolmogorov-Smirnov distribution for a sample size of N.
    
    properties(SetAccess = protected)
        N     % Sample size
    end

    properties(Constant)    
        minN = 1;
        sJarFileName = 'KolmoJava1_7Class.jar';
    end

    methods
        
        function obj=KolmSmir(varargin)   % Constructor
            obj=obj@dContinuous('KolmSmir');
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            obj.NDistParms = 1;
            p = mfilename('fullpath');
            mypath = fileparts(p);
            javaaddpath([mypath '\' obj.sJarFileName]);
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('KolmSmir:Constructor', ...
                        'Too many arguments passed to KolmSmir constructor.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.N = VerifyIntegerGE(obj,obj.minN,newparmvalues(1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            if ~(ParmCodes(1) == 'f')
                obj.ResetParms(obj.N+1);
            end
        end
        
        function []=ReInit(obj)
            assert((obj.N>=obj.minN)&&(iswholenumber(obj.N)),'KolmSmir N must be an integer >= 2.');
            obj.LowerBound = 0;
            obj.UpperBound = 1;
            obj.Initialized = true;  % Needed so that CDF can be called
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.GT2Real(obj.minN,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(obj.minN,Reals(1));
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for i=1:numel(X)
                if InBounds(i)
                   thiscdf(i) =  KolmogorovSmirnovDist.cdf(obj.N,X(i));
                end
            end
        end
        
    end  % methods
    
end  % class KolmSmir

