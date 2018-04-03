classdef MonotoneTrans < dTransOf1  % NWJEFF: Not documented.
    % MonotoneTrans(BasisRV,TransFunc,InverseFunc) produces any monotonic transformation of the original BasisRV.
    % TransFunc is the function that does the transformation
    % InverseFunc is the function that does the inverse of the transformation
    % Note that TransFunc and InverseFunc are fixed--i.e., they have no adjustable parameters.
    %
    % Here are a few examples:
    %
    % f = @(x) x+5
    % finv=@(x) x-5
    % md = MonotoneTrans(Uniform(0,1),f,finv)
    %
    % f2=@(x) sqrt(x)
    % f2inv=@(x) x.^2
    % md2=MonotoneTrans(Uniform(10,12),f2,f2inv)
    %
    % fRolf=@(x) x./(1-x)
    % fInvRolf=@(x) x./(1+x)
    % distRolf = MonotoneTrans(Uniform(.1,.9),fRolf,fInvRolf);
    
    properties(SetAccess = protected)
        TransFunc, InverseFunc, Increasing
    end
    
    methods
        
        function obj=MonotoneTrans(varargin)
            obj=obj@dTransOf1('MonotoneTrans');
            obj.NTransParms = 0;
            obj.TransParmCodes = '';
            switch nargin
                case 0
                case 3
                    BuildMyBasis(obj,varargin{1});
                    obj.TransFunc = varargin{2};
                    obj.InverseFunc = varargin{3};
                    obj.DistType = obj.BasisRV.DistType;
                    % NWJEFF: No parameters allowed.
                    obj.NDistParms = obj.BasisRV.NDistParms;
                    obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
                    ResetParms(obj,[obj.BasisRV.ParmValues]);
                otherwise
                    ME = MException('MonotoneTrans:Constructor', ...
                        'MonotoneTrans constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = ['MonotoneTrans(' obj.BasisRV.StringName ',' func2str(obj.TransFunc)  ',' func2str(obj.InverseFunc) ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransOf1(obj,newparmvalues);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            obj.ResetParms(obj.BasisRV.ParmValues);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.Increasing = obj.TransFunc(obj.BasisRV.UpperBound) > obj.TransFunc(obj.BasisRV.LowerBound);
            if obj.Increasing
                obj.LowerBound = obj.TransFunc(obj.BasisRV.LowerBound);
                obj.UpperBound = obj.TransFunc(obj.BasisRV.UpperBound);
            else
                obj.LowerBound = obj.TransFunc(obj.BasisRV.UpperBound);
                obj.UpperBound = obj.TransFunc(obj.BasisRV.LowerBound);
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = obj.TransFunc(PreTransX);
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = obj.InverseFunc(TransX);
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = [];
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = [];
        end
        
        function thisval = nIthValue(obj,Ith)
            thisval = obj.TransFunc(obj.BasisRV.nIthValue(Ith));
        end
        
        function thispdf=PDF(obj,X)
            thispdf = PDF@dContinuous(obj,X); % dEither(obj,X);
        end
        
        function thisval = PDFScaleFactor(obj,~)
            thisval = 1;  % Not needed because PDF is obtained by numerical differentiation of CDF.
        end
        
    end  % methods
    
end  % class MonotoneTrans


