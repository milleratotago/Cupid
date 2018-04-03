classdef (Abstract) dTransOf1 < dContinuous % dEither
    % An abstract class to encapsulate code common across many monotonic transformations
    % of a single RV, including AddTrans, MultTrans, LinearTrans, PowerTrans, ExpTrans, LogTrans

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
        NTransParms      % Number of parameters used by the transformation, e.g. 1 for AddTrans, 0 for LogTrans
        TransParmCodes   % ParmCodes string for the transformation parameters (usually 'r')
        TransReverses    % Set to true if the transformation reverses the mapping,
                         % i.e. if the smallest basis values are mapped to
                         % the largest values in the transformed distribution.
                         % For example, multiplying by a negative number
    end
    
    methods(Abstract)
        parmvals = TransParmValues(obj)
        TransX = PreTransToTrans(obj,PreTransX)
        PreTransX = TransToPreTrans(obj,TransX)
        TransReals = TransParmsToReals(obj,Parms,~)
        TransParms = TransRealsToParms(obj,Reals,~)
        thisval = PDFScaleFactor(obj,X)
    end
    
    methods
        
        function obj=dTransOf1(FamName)
            obj@dContinuous(FamName); % dEither(FamName);
            obj.TransReverses = false;
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues)
            obj.Initialized = false;
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NValues = obj.BasisRV.NValues;
            obj.NDistParms = obj.BasisRV.NDistParms + obj.NTransParms;
            obj.DefaultParmCodes = [obj.BasisRV.DefaultParmCodes obj.TransParmCodes];
            % ReInit(obj);
        end
        
        function BuildMyBasis(obj,s)
            obj.BasisRV = s;
            assert(obj.BasisRV.Initialized,['Error initializing ' obj.ThisFamilyName ' BasisRV.']);
            obj.DistType = obj.BasisRV.DistType;
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.ThisFamilyName '(' obj.BasisRV.StringName];
            parms = ParmValues(obj);
            for iParm=1:obj.NTransParms
                obj.StringName = [obj.StringName ',' num2str(parms(obj.BasisRV.NDistParms+iParm))];
            end
            obj.StringName = [obj.StringName ')'];
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.BasisRV.ParmsToReals(Parms(1:obj.BasisRV.NDistParms)) TransParmsToReals(obj,Parms)];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.BasisRV.RealsToParms(Reals(1:obj.BasisRV.NDistParms)) TransRealsToParms(obj,Reals)];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.BasisRV.ParmValues TransParmValues(obj)];
        end
        
        function x = XsToPlot(obj)
            x = PreTransToTrans(obj,XsToPlot(obj.BasisRV));
        end
        
        function thisval=NearestLegal(obj,X)
            PreTrans = TransToPreTrans(obj,X);
            PreTrans = obj.BasisRV.NearestLegal(PreTrans);
            thisval = PreTransToTrans(obj,PreTrans);
        end
        
        function thisval=LegalValue(obj,X)
            PreTrans = TransToPreTrans(X);
            thisval = obj.BasisRV.LegalValue(PreTrans);
        end
        
        function thisval=IthValue(obj,I)
            thisval = PreTransToTrans(obj,obj.BasisRV.IthValue(I));
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            PreTransX = TransToPreTrans(obj,X(InBounds));
            thispdf(InBounds) = PDF(obj.BasisRV,PreTransX) .* PDFScaleFactor(obj,X(InBounds));
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            PreTransX = TransToPreTrans(obj,X(InBounds));
            if obj.TransReverses
                thiscdf(InBounds) = 1 - CDF(obj.BasisRV,PreTransX);
            else
                thiscdf(InBounds) = CDF(obj.BasisRV,PreTransX);
            end
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            if obj.TransReverses
                P = 1 - P;
            end
            PreTransX = InverseCDF(obj.BasisRV,P);
            thisval = PreTransToTrans(obj,PreTransX);
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = PreTransToTrans(obj,Random(obj.BasisRV,varargin{:}));
        end
        
    end  % methods
    
end  % class dTransOf1


