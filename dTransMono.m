classdef (Abstract) dTransMono < dEither
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
        PDFScaleFactorKnown
    end
    
    methods(Abstract)
        TransX = PreTransToTrans(obj,PreTransX)
        PreTransX = TransToPreTrans(obj,TransX)
    end
    
    methods
        
        function obj=dTransMono(FamName,BasisDist)
            obj@dEither(FamName);
            obj.BasisRV = BasisDist;
            obj.DistType = obj.BasisRV.DistType;
            obj.NTransParms = 0;
            obj.NDistParms = obj.BasisRV.NDistParms;
            obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
        end

        function []=AddParms(obj,NTransParms,TransParmCodes)
            obj.NTransParms = NTransParms;
            obj.DefaultParmCodes = [obj.DefaultParmCodes TransParmCodes];
            obj.NDistParms = obj.NDistParms + NTransParms;
        end
        
        function PerturbParms(obj,ParmCodes)    % Default for descendants without parms
            obj.BasisRV.PerturbParms(ParmCodes);
            obj.ResetParms(obj.BasisRV.ParmValues);
        end
        
        function parmvals = TransParmValues(obj)    % Default for descendants without parms
            parmvals = [];
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)    % Default for descendants without parms
            TransReals = [];
        end
        
        function Parms=TransRealsToParms(obj,Reals)    % Default for descendants without parms
            Parms = [];
        end
        
        function []=ResetParms(obj,newparmvalues)    % Default for descendants without parms
            ClearBeforeResetParms(obj)
            obj.Initialized = true;
            obj.BasisRV.ResetParms(newparmvalues(1:obj.BasisRV.NDistParms));
            obj.NValues = obj.BasisRV.NValues;
            % obj.TransReverses = obj.PreTransToTrans(obj.BasisRV.UpperBound) < obj.PreTransToTrans(obj.BasisRV.LowerBound);
            % ReInit(obj);  % Descendants with parms must set them before ReInit
        end
        
         function []=ReInit(obj)
            obj.Initialized = true;
            switch obj.DistType
                case 'd'
                    obj.MakeTables;
                case 'c'
                    obj.InitContin;
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
       function []=InitContin(obj)
            if obj.TransReverses
                obj.LowerBound = obj.PreTransToTrans(obj.BasisRV.UpperBound);
                obj.UpperBound = obj.PreTransToTrans(obj.BasisRV.LowerBound);
            else
                obj.LowerBound = obj.PreTransToTrans(obj.BasisRV.LowerBound);
                obj.UpperBound = obj.PreTransToTrans(obj.BasisRV.UpperBound);
            end
            % Improve bounds to avoid numerical errors
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
        end
        
        function BuildMyBasis(obj,s)
            obj.BasisRV = s;
            assert(obj.BasisRV.Initialized,['Error initializing ' obj.FamilyName ' BasisRV.']);
            obj.DistType = obj.BasisRV.DistType;
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' obj.BasisRV.StringName];
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
        
        function MakeTables(obj)
            obj.NValues = obj.BasisRV.NValues;
            obj.DiscreteX = PreTransToTrans(obj,obj.BasisRV.DiscreteX);
            obj.DiscretePDF = obj.BasisRV.DiscretePDF;
            obj.DiscreteCDF = obj.BasisRV.DiscreteCDF;
            if obj.TransReverses
                obj.DiscreteX = flip(obj.DiscreteX);
                obj.DiscretePDF = flip(obj.DiscretePDF);
                % 1-CDF reflects the probability of scores previously > than the current value,
                % and PDF adds in the probability of scores == the current value.
                obj.DiscreteCDF = 1 - flip(obj.DiscreteCDF) + obj.DiscretePDF;
            end
            obj.SetBinEdges;
            obj.LowerBound = obj.DiscreteXmin(1);
            obj.UpperBound = obj.DiscreteXmax(end);
        end
        
        function thispdf=PDF(obj,X)
            if obj.DistType=='d'
                thispdf = PDF@dDiscrete(obj,X);
                return;
            end
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if (obj.DistType=='c') && obj.PDFScaleFactorKnown
                PreTransX = TransToPreTrans(obj,X(InBounds));
                thispdf(InBounds) = PDF(obj.BasisRV,PreTransX) .* PDFScaleFactor(obj,X(InBounds));
            elseif obj.DistType=='c'
                thispdf(InBounds) = PDF@dContinuous(obj,X(InBounds));
            else
                ME = MException('dTransMono:PDF','Unrecognized distribution type.');
                throw(ME);
            end
        end
        
        function thiscdf=CDF(obj,X)
            if obj.DistType=='d'
                thiscdf = CDF@dDiscrete(obj,X);
                return;
            end
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
            thiscdf(thiscdf<0) = 0;
            thiscdf(thiscdf>1) = 1;
        end
        
        %         function thisval=InverseCDF(obj,P)
        %             assert(obj.Initialized,UninitializedError(obj));
        %             assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
        %             if obj.TransReverses
        %                 P = 1 - P;
        %             end
        %             PreTransX = InverseCDF(obj.BasisRV,P);
        %             thisval = PreTransToTrans(obj,PreTransX);
        %         end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            %             if obj.DistType=='d'
            %                 thisval = Random@dDiscrete(obj,varargin{:});
            %             else
            thisval = PreTransToTrans(obj,Random(obj.BasisRV,varargin{:}));
            %             end
        end
        
        function thisval = PDFScaleFactor(obj,X)  % Default for unknown
            thisval = [];
        end

    end  % methods
    
end  % class dTransMono


