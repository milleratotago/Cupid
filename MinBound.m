classdef MinBound < dContinuous % dEither
    % MinBound(BasisRV1,BasisRV2) creates a random variable corresponding to the minimum of two
    % nonindependent racer distributions based on the function
    %          Fm(t) = F1(t) + F2(t)
    % where F1 and F2 are the finishing times of the two racers and Fm is this minimum distribution.
    
    properties(SetAccess = protected)
        BasisRV1, BasisRV2
    end
    
    methods
        
        function obj=MinBound(varargin)
            obj=obj@dContinuous('MinBound');
            switch nargin
                case 0
                case 2
                    obj.BasisRV1 = varargin{1};
                    obj.BasisRV2 = varargin{2};
                    assert(obj.BasisRV1.Initialized&&obj.BasisRV2.Initialized,['Error initializing ' obj.FamilyName ' BasisRV1/RV2.']);
                    if obj.BasisRV1.DistType == obj.BasisRV2.DistType
                        obj.DistType = obj.BasisRV1.DistType;
                    else
                        warning('MinBound racers must be either both continuous or both discrete---not one of each.');
                    end
                    obj.DistType = obj.BasisRV1.DistType;
                    if (obj.BasisRV1.DistType=='c') && (obj.BasisRV2.DistType=='c')
                        obj.DistType = 'c';
                    else
                        assert(false,'MinBound can only handle continuous Basis distributions (so far)');
                    end
                    obj.NDistParms = obj.BasisRV1.NDistParms + obj.BasisRV2.NDistParms;
                    obj.DefaultParmCodes = [obj.BasisRV1.DefaultParmCodes obj.BasisRV2.DefaultParmCodes];
                    ResetParms(obj,[obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
                otherwise
                    ME = MException('MinBound:Constructor', ...
                        'MinBound constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' obj.BasisRV1.StringName ',' obj.BasisRV2.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj)
            obj.Initialized = false;
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
            obj.LowerBound = min([obj.BasisRV1.LowerBound obj.BasisRV2.LowerBound]);
            obj.UpperBound = min([obj.BasisRV1.UpperBound obj.BasisRV2.UpperBound]);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.BasisRV1.ParmsToReals(Parms(1:obj.BasisRV1.NDistParms)) obj.BasisRV2.ParmsToReals(Parms(obj.BasisRV1.NDistParms+1:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.BasisRV1.RealsToParms(Reals(1:obj.BasisRV1.NDistParms)) obj.BasisRV2.RealsToParms(Reals(obj.BasisRV1.NDistParms+1:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues];
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            if obj.DistType == 'd'
                thiscdf = LRV.CDF(X);
            else
                FX1 = obj.BasisRV1.CDF(X);
                FX2 = obj.BasisRV2.CDF(X);
                thiscdf = FX1 + FX2;
                thiscdf(thiscdf>1) = 1;
            end
        end
        
    end  % methods
    
end  % class MinBound

% function thisval=MinBoundRV.InitDiscrete;
%        obj.DistType = 'd';
% MaxN = obj.BasisRV1.NValues + obj.BasisRV2.NValues;
% %Writeln('MaxN=',MaxN);
% LRV.Init(MaxN);
% Ptr1 = 1;
% Ptr2 = 1;
% X1 = obj.BasisRV1.ithValue(1);
% X2 = obj.BasisRV2.ithValue(1);
% iList = 0;
% Lastthiscdf = 0;
% Repeat
%    if X1 <= X2 Then ThisX = X1 else ThisX = X2;
%    if X1 <= ThisX
%       Inc(Ptr1);
%       X1 = obj.BasisRV1.ithValue(Ptr1);
%       end
%    if X2 <= ThisX
%       Inc(Ptr2);
%       X2 = obj.BasisRV2.ithValue(Ptr2);
%       end
%    Thisthiscdf = obj.BasisRV1.CDF(ThisX) + obj.BasisRV2.CDF(ThisX);
%    if ThisCDF > 1 Then Thisthiscdf = 1;
%    Inc(iList);
%    LRV.XVal(iList) = ThisX;
%    LRV.PDFVal(iList) = ThisCDF - LastCDF;
%    Lastthiscdf = ThisCDF;
% %   Writeln('Position ',iList,' storing ',ThisX:6:2,' with PDF=',LRV.PDFVal(iList):5:3);
% %   Writeln(X1:5:1,' ',X2:5:1,' ',ThisCDF:7:4,LastCDF:7:4);
% %   Readln;
%  Until (ThisCDF >= obj.CDFNearlyOne) or ((Ptr1 = obj.BasisRV1.NValues) and (Ptr2 = obj.BasisRV2.NValues));
% LRV.Summarize(iList);
% %Writeln('Passed summarize.');
% %With LRV Do Writeln('NValues=',NValues,', obj.LowerBound=',obj.LowerBound:8:2,', obj.UpperBound=',obj.UpperBound:8:2);
% NValues = LRV.NValues;
% obj.LowerBound = LRV.obj.LowerBound;
% obj.UpperBound = LRV.obj.UpperBound;
%   obj.Initialized = true;
% end
