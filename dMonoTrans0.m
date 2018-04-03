classdef dMonoTrans0 < dContinuous  % dEither
    % This is a base class to handle parameter-free, invertable, monotonic transformations of a BasisRV.
    % It may be better to use dTransOf1 if the transformation is differentiable so that PDFScaleFactor can be implemented.
    
    % Notes:
    %   PDF will often fail with discrete RVs because of floating point imprecision.

% Copyright (C) 2018 Jeffrey Owen Miller
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program (00License.txt). If not, see 
%     <http://www.gnu.org/licenses/>.
%

    properties(SetAccess = protected)
        BasisRV
        TransReverses    % Set to true if the transformation reverses the mapping,
                         % i.e. if the smallest basis values are mapped to
                         % the largest values in the transformed distribution.
    end
    
    methods(Abstract)
        
        Trans = PreTransToTrans(obj,PreTrans)
        PreTrans = TransToPreTrans(obj,Trans)
        
    end  % Abstract methods

    methods
        
        function obj=dMonoTrans0(FamName,BasisDist)
            obj=obj@dContinuous(FamName);
            obj.BasisRV = BasisDist;
            obj.NDistParms = obj.BasisRV.NDistParms;
            obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.ThisFamilyName '(' obj.BasisRV.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.BasisRV.CheckBeforeResetParms(newparmvalues);
            obj.BasisRV.ResetParms(newparmvalues);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes);
            obj.ResetParms(obj.BasisRV.ParmValues);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.NValues = obj.BasisRV.NValues;
            obj.LowerBound = PreTransToTrans(obj,obj.BasisRV.LowerBound);
            obj.UpperBound = PreTransToTrans(obj,obj.BasisRV.UpperBound);
            % Improve bounds to avoid numerical errors
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
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
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if obj.DistType == 'd'
                PreTransX = TransToPreTrans(obj,X(InBounds));
                thispdf(InBounds) = PDF(obj.BasisRV,PreTransX);
            else
                thispdf(InBounds) = PDF@dContinuous(obj,X(InBounds));
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            PreTransX = TransToPreTrans(obj,X(InBounds));
            thiscdf(InBounds) = CDF(obj.BasisRV,PreTransX);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            PreTransX = InverseCDF(obj.BasisRV,P);
            thisval = PreTransToTrans(obj,PreTransX);
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = PreTransToTrans(obj,Random(obj.BasisRV,varargin{:}));
        end
        
        function thisval=LegalValue(obj,X)
            temp = obj.TransToPreTrans(obj,X);
            thisval = obj.BasisRV.LegalValue(temp);
        end
        
        function thisval=NearestLegal(obj,X)
            temp = obj.TransToPreTrans(X);
            temp = obj.BasisRV.NearestLegal(temp);
            thisval = obj.PreTransToTrans(temp);
        end
        
        function thisval = nIthValue(obj,Ith)
            thisval = PreTransToTrans(obj.BasisRV.nIthValue(Ith));
        end
        
    end  % methods
    
end  % class MonoTrans2

