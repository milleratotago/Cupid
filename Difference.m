classdef Difference < dContinuous  % dEither    % Only continuous so far (old discrete saved at end).  Should handle differences of discrete RVs via ListRV mechanism.
    % Difference(BasisRV1,BasisRV2) creates a random variable that is
    %  the difference between two independent basis random variables,
    %  BasisRV1-BasisRV2
    
    properties(SetAccess = protected)
        BasisRV1, BasisRV2
    end
    
    methods
        
        function obj=Difference(varargin)
            obj=obj@dContinuous('Difference');
            switch nargin
                case 0
                case 2
                    obj.BasisRV1 = varargin{1};
                    obj.BasisRV2 = varargin{2};
                    if (obj.BasisRV1.DistType=='c') && (obj.BasisRV2.DistType=='c')
                        obj.DistType = 'c';
                    else
                        assert(false,'Difference can only handle continuous Basis distributions (so far)');
                    end
                    obj.NDistParms = obj.BasisRV1.NDistParms + obj.BasisRV2.NDistParms;
                    obj.DefaultParmCodes = [obj.BasisRV1.DefaultParmCodes obj.BasisRV2.DefaultParmCodes];
                    ResetParms(obj,[obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
                otherwise
                    ME = MException('Difference:Constructor', ...
                        'Difference constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            assert(obj.Initialized,UninitializedError(obj));
            obj.StringName = ['Difference(' obj.BasisRV1.StringName ',' obj.BasisRV2.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.BasisRV1.ResetParms(newparmvalues(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.ResetParms(newparmvalues(obj.BasisRV1.NDistParms+1:end));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV1.PerturbParms(ParmCodes(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.PerturbParms(ParmCodes(obj.BasisRV1.NDistParms+1:end));
            obj.ResetParms([obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = obj.BasisRV1.LowerBound - obj.BasisRV2.UpperBound;
            obj.UpperBound = obj.BasisRV1.UpperBound - obj.BasisRV2.LowerBound;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    thispdf(iel)=integral(@(x) obj.BasisRV1.PDF(X(iel)+x).*obj.BasisRV2.PDF(x),obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound);
                end
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            %            assert(obj.Initialized,UninitializedError(obj));
            Reals = [obj.BasisRV1.ParmsToReals(Parms(1:obj.BasisRV1.NDistParms)) obj.BasisRV2.ParmsToReals(Parms(obj.BasisRV1.NDistParms+1:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            %            assert(obj.Initialized,UninitializedError(obj));
            Parms = [obj.BasisRV1.RealsToParms(Reals(1:obj.BasisRV1.NDistParms)) obj.BasisRV2.RealsToParms(Reals(obj.BasisRV1.NDistParms+1:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues];
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Mean - obj.BasisRV2.Mean;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Variance + obj.BasisRV2.Variance;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Random(varargin{:}) - obj.BasisRV2.Random(varargin{:});
        end
        
    end  % methods
    
end  % class Difference


% **************** ddiffer.pas
%
% DisDifferenceRV = Class(DifferenceRV) % This is the discrete 1.
%    LRV : ListRV;
%    Constructor Create;
%    function thisval=Init;    % All the BasisRV distributions must already have been
%                        initialized before either init is called.
%    function thisval=PDF(X);
%    function thisval=CDF(X);
%    function thisval=NearestLegal(X);
%    function thisval=IthValue(I );
%    end
%
% %*********Discrete Difference Random Variables:****************************}
%
%
% % An almost identical routine is used to set up discrete Convolution RVs.
% function thisval=DisDifferenceRV.Init;
% Var iVal1, iVal2, iValY, N ;
%     X1, X2, Y, Pr1, Pr2;
% Inherited Init;
% if (BasisRV1.DistType <> Discrete) or (BasisRV2.DistType <> Discrete)
%    warning('DisDifferenceRV.Init should only be called with discrete BasisRV distributions.');
%    Exit;
%    end
% N = BasisRV1.NValues*BasisRV2.NValues;
% LRV.Init(N);
% iValY = 0;
% for iVal1 = 1:BasisRV1.NValues Do Begin
%    X1 = BasisRV1.ithValue(iVal1);
%    Pr1 = BasisRV1.PDF(X1);
%    for iVal2 = 1:BasisRV2.NValues Do Begin
%       X2 = BasisRV2.ithValue(iVal2);
%       Pr2 = BasisRV2.PDF(X2);
%       Y = X1 - X2;  % Here is the Difference step.
%       Inc(iValY);
%       LRV.XVal(iValY) = Y;
%       LRV.PDFVal(iValY) = Pr1*Pr2;
%       end
%    end
% LRV.Summarize(N);
% NValues = LRV.NValues;
% Initialized = true;
% UCDResult = OK;
% DistType = Discrete;
% end
%
% function thisval=DisDifferenceRV.PDF(X);
% thispdf = LRV.PDF(X);
% end
%
% function thisval=DisDifferenceRV.CDF(X);
% thiscdf = LRV.CDF(X);
% end
%
% function thisval=DisDifferenceRV.IthValue(I );
% IthValue = LRV.IthValue(I);
% end
%
% function thisval=DisDifferenceRV.NearestLegal(X);
% NearestLegal = LRV.NearestLegal(X);
% end
