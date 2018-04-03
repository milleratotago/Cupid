classdef OrderIID < dContinuous % dEither
    % OrderIID(Order,SampleSize,BasisRV): The Order'th order statistic of SampleSize=N IID RVs from the BasisRV.
    
    % Notes:
    % o By default Order & SampleSize are NOT adjusted when parameter fitting.
    % o CDF is computed directly but PDF is obtained from CDF.
    % o Has no ReInit because ???
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        Order, SampleSize, BasisRV, CDFLoop1, CDFLoop2, BigOrder
    end
    
    methods
        
        function obj=OrderIID(varargin)   % Constructor
            obj=obj@dContinuous('OrderIID'); % dEither('OrderIID');
            switch nargin
                case 0
                case 3
                    obj.BasisRV = varargin{3};
                    assert(obj.BasisRV.Initialized,['Error initializing ' obj.ThisFamilyName ' as BasisRV for OrderIID.']);
                    obj.DistType = obj.BasisRV.DistType;
                    obj.NDistParms = 2 + obj.BasisRV.NDistParms;
                    obj.DefaultParmCodes = ['ff' obj.BasisRV.DefaultParmCodes];
                    ResetParms(obj,[varargin{1:2} obj.BasisRV.ParmValues]);
                otherwise
                    ME = MException('OrderIID:Constructor', ...
                        'Illegal number of arguments passed to OrderIID constructor.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = ['OrderIID(' num2str(obj.Order) ',' num2str(obj.SampleSize) ',' obj.BasisRV.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Initialized = false;
            CheckBeforeResetParms(obj,newparmvalues);
            obj.Order = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.SampleSize = VerifyIntegerGE(obj,obj.Order,newparmvalues(2));
            obj.BigOrder = obj.Order >= floor( (obj.SampleSize + 1) / 2 );
            if obj.BigOrder
                obj.CDFLoop1 = obj.Order;
                obj.CDFLoop2 = obj.SampleSize;
            else
                obj.CDFLoop1 = 0;
                obj.CDFLoop2 = obj.Order - 1;
            end
            
            obj.BasisRV.ResetParms(newparmvalues(3:end));
            obj.NValues = obj.BasisRV.NValues;
            % Finding bounds is difficult when the basis distribution is unbounded
            % (e.g., OrderIID(1,10,Normal(0,1)) ).  Some cluges used in Pascal version.
            %      Adjust = BasisRV.SD / 2;
            %      obj.LowerBound = BasisRV.obj.LowerBound;
            %      obj.UpperBound = BasisRV.obj.UpperBound;
            %      While CDF(obj.LowerBound) > obj.CDFNearlyZero Do obj.LowerBound = obj.LowerBound - Adjust;
            %      While CDF(obj.UpperBound) < obj.CDFNearlyOne Do obj.UpperBound = obj.UpperBound + Adjust;
            obj.LowerBound = obj.BasisRV.LowerBound;
            obj.UpperBound = obj.BasisRV.UpperBound;
            obj.Initialized = true;
            if obj.NameBuilding
                BuildMyName(obj);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV.PerturbParms(ParmCodes(3:end));
            obj.ResetParms(obj.ParmValues);
        end
        
        function parmvals = ParmValues(obj,varargin)
            parmvals = [obj.Order obj.SampleSize obj.BasisRV.ParmValues];
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % Improve here and next procedure by constraining Order <= SampleSize.
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.GT2Real(1,Parms(2)) obj.BasisRV.ParmsToReals(Parms(3:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2GT(1,Reals(2)) obj.BasisRV.RealsToParms(Reals(3:end))];
        end
        
        function thisval=LegalValue(obj,X)
            thisval = LegalValue(obj.BasisRV,X);
        end
        
        function thisval=NearestLegal(obj,X)
            thisval = NearestLegal(obj.BasisRV,X);
        end
        
        function thisval=nIthValue(obj,I)
            thisval = nIthValue(obj.BasisRV,I);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    Fx = CDF(obj.BasisRV,X(iel));
                    if Fx <= 0
                        thiscdf(iel) = 0;
                    elseif Fx >= 1
                        thiscdf(iel) = 1;
                    else
                        LnFx = log(Fx);
                        Ln1_Fx = log(1 - Fx);
                        FxSum = 0;
                        N = nchoosek(obj.SampleSize,obj.CDFLoop1);
                        for Idx = obj.CDFLoop1:obj.CDFLoop2
                            One_FxTerm = Ln1_Fx * (obj.SampleSize - Idx);
                            FxTerm = LnFx * Idx;
                           FxSum = exp(FxTerm + One_FxTerm) * N + FxSum;
                            N = N * (obj.SampleSize - Idx) / (Idx + 1);
                        end
                        if obj.BigOrder
                            thiscdf(iel) = FxSum;
                        else
                            thiscdf(iel) = 1 - FxSum;
                        end
                    end
                end
            end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval=zeros(varargin{:});
            for i=1:numel(thisval)
                X = Random(obj.BasisRV,obj.SampleSize,1);
                X = sort(X);
                thisval(i) = X(obj.Order);
            end
        end
        
    end  % methods
    
end  % class OrderIID

