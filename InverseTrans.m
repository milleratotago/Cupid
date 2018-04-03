classdef InverseTrans < dTransOf1  % NWJEFF: change to dMonoTrans1
    % InverseTrans(BasisRV): Inverse (1/X) transformation of BasisRV, which must be either all positive or all negative values.
    
    properties(SetAccess = protected)
    end
    
    methods
        
        function obj=InverseTrans(varargin)
            obj=obj@dTransOf1('InverseTrans');
            obj.NTransParms = 0;
            obj.TransParmCodes = '';
            switch nargin
                case 0
                case 1
                    BuildMyBasis(obj,varargin{1});
                    obj.DistType = obj.BasisRV.DistType;
                    obj.NDistParms = obj.BasisRV.NDistParms;
                    obj.DefaultParmCodes = obj.BasisRV.DefaultParmCodes;
                    ResetParms(obj,obj.BasisRV.ParmValues);
                otherwise
                    ME = MException('InverseTrans:Constructor', ...
                        'InverseTrans constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
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
            assert(PDF(obj.BasisRV,0)==0,'InverseTrans BasisRV must have PDF(0)=0.');
            assert(obj.BasisRV.LowerBound*obj.BasisRV.UpperBound>0,'InverseTrans BasisRV values must be all positive or all negative.');
            obj.Initialized = true;
            obj.LowerBound = 1/obj.BasisRV.LowerBound;
            obj.UpperBound = 1/obj.BasisRV.UpperBound;
            if obj.BasisRV.LowerBound >= 0   % obj.BasisRV is all positive.
                obj.TransReverses = true;
                obj.LowerBound = 1 / obj.BasisRV.UpperBound;
                obj.UpperBound = 1 / obj.BasisRV.LowerBound;
            elseif obj.BasisRV.UpperBound <= 0   % obj.BasisRV is all negative.
                obj.TransReverses = false;
                obj.LowerBound = 1 / obj.BasisRV.UpperBound;
                obj.UpperBound = 1 / obj.BasisRV.LowerBound;
            else    % obj.BasisRV has some positive values & some negative.
                % It is probably possible to handle this case but I have not done so.
                ME = MException('InverseTrans:Init', ...
                    'InverseTrans Basis RV cannot include both positive and negative values.');
                throw(ME);
            end
            % Improve bounds to avoid numerical errors
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = TransParmValues(obj)
            parmvals = [];
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = ones(size(PreTransX)) ./ PreTransX;
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = ones(size(TransX)) ./ TransX;
        end
        
        function TransReals = TransParmsToReals(obj,Parms,~)
            TransReals = [];
        end
        
        function TransParms = TransRealsToParms(obj,Reals,~)
            TransParms = [];
        end
        
        function thisval = PDFScaleFactor(obj,X)
            thisval=ones(size(X));
            switch obj.DistType
                case 'c'
                    PreTransX = ones(size(X)) ./ X;
                    thisval = thisval .* PreTransX.^2;
            end
        end
        
        function thisval = nIthValue(obj,Ith)
            thisval = TransX(obj.BasisRV.nIthValue(Ith));
        end
        
    end  % methods
    
end  % class InverseTrans


