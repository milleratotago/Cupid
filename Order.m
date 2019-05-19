classdef Order < dEither
    % Order(Order,BasisRV1,BasisRV2,BasisRV3,...): The Order'th order statistic of the independent but non-identical BasisRVs.
    
    % Notes:
    % o By default Order is NOT adjusted when parameter fitting.
    % o CDF is computed directly but PDF is obtained from CDF.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        OrderK, SampleSize, BasisRV, CDFLoop1, CDFLoop2, BigOrder,
        CumNParms, DistParmP, DistParmCodes, PreviousParmCodes
    end
    
    methods
        
        function obj=Order(varargin)   % Constructor
            obj=obj@dEither('Order');
            switch nargin
                case 0
                case {1, 2}
                    ME = MException('Order:Constructor', ...
                        'Order constructor needs 0 or 3+ arguments.');
                    throw(ME);
                otherwise
                    Setup(obj,varargin(:));
                    ResetParms(obj,[obj.ParmValues]);
            end
        end
        
        function Setup(obj,s)
            obj.SampleSize = numel(s) - 1;
            obj.Initialized = false;
            obj.BasisRV = cell(obj.SampleSize,1);
            obj.CumNParms = zeros(obj.SampleSize,1);
            obj.DistParmP = char(obj.SampleSize);
            obj.DistParmCodes = cell(obj.SampleSize,1);
            obj.PreviousParmCodes = '-1';   % Reset to an impossible value so that new ones will be computed by CheckNewParmCodes.
            obj.OrderK = s{1};
            for iDist=1:obj.SampleSize
                nextptr = iDist + 1;
                obj.BasisRV{iDist} = s{nextptr};
                if ~obj.BasisRV{iDist}.Initialized
                    error(['Unable to initialize Order BasisRV number ' num2str(iDist)]);
                end
            end
            
            % Determine distribution type:
            OneIsDiscrete = false;
            OneIsContinuous = false;
            OneIsMixed = false;
            for I = 1:obj.SampleSize
                if obj.BasisRV{I}.DistType == 'c'; OneIsContinuous = true; end
                if obj.BasisRV{I}.DistType == 'd'; OneIsDiscrete = true; end
                if obj.BasisRV{I}.DistType == 'm'; OneIsMixed = true; end
            end
            if OneIsContinuous && ~(OneIsDiscrete || OneIsMixed)
                obj.DistType = 'c';
            elseif OneIsDiscrete && ~(OneIsContinuous || OneIsMixed)
                obj.DistType = 'd';
            else
                error('Order can only handle all-continuous or all-discrete Basis distributions (so far)');
            end
            
            % Assemble DefaultParmCodes & count NDistParms & CumNParms
            obj.NDistParms = 1;       % Order
            obj.DefaultParmCodes = 'f';
            for iDist=1:obj.SampleSize
                obj.DefaultParmCodes = [obj.DefaultParmCodes obj.BasisRV{iDist}.DefaultParmCodes];
                obj.NDistParms = obj.NDistParms + obj.BasisRV{iDist}.NDistParms;
                obj.CumNParms(iDist) = obj.NDistParms;
            end
            
        end
        
        function BuildMyName(obj)
            sTemp = ['Order(' num2str(obj.OrderK)];
            for iDist=1:obj.SampleSize
                sTemp = [sTemp ',' obj.BasisRV{iDist}.StringName];
            end
            obj.StringName = [sTemp ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj);
            obj.Initialized = false;
            obj.OrderK = VerifyIntegerInRange(obj,1,obj.SampleSize,newparmvalues(1));
            obj.BigOrder = obj.OrderK >= floor( (obj.SampleSize + 1) / 2 );
            if obj.BigOrder
                obj.CDFLoop1 = obj.OrderK;
                obj.CDFLoop2 = obj.SampleSize;
            else
                obj.CDFLoop1 = 0;
                obj.CDFLoop2 = obj.OrderK - 1;
            end
            
            ParmsUsed = 1;  % OrderK
            for iDist = 1:obj.SampleSize
                obj.BasisRV{iDist}.ResetParms(newparmvalues(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms));
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
            
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            ParmsUsed = 1;  % OrderK
            for iDist = 1:obj.SampleSize
                ThisDistCodes = ParmCodes(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms);
                obj.BasisRV{iDist}.PerturbParms(ThisDistCodes);
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
            obj.ResetParms(obj.ParmValues);
        end
        
        function ReInit(obj)
            % Assume all distributions have been initialized
            
            Extremes = zeros(obj.SampleSize,1);
            for iDist=1:obj.SampleSize
                Extremes(iDist) = obj.BasisRV{iDist}.LowerBound;
            end
            Extremes = sort(Extremes);
            obj.LowerBound = Extremes(obj.OrderK);
            for iDist=1:obj.SampleSize
                Extremes(iDist) = obj.BasisRV{iDist}.UpperBound;
            end
            Extremes = sort(Extremes);
            obj.UpperBound = Extremes(obj.OrderK);
            obj.Initialized = true;
            if obj.DistType == 'd'
                obj.MakeTables;
            end
            % obj.NValues = obj.BasisRV.NValues;  % dEither
            if obj.NameBuilding
                BuildMyName(obj);
            end
            
        end
        
        function MakeTables(obj)
            Xs = [];
            for iDist=1:obj.SampleSize
                Xs = [Xs obj.BasisRV{iDist}.DiscreteX];
            end
            obj.DiscreteX = unique(Xs);
            obj.DistType = 'n';  % Temporarily turn off 'd' so we can use CDF
            obj.DiscreteCDF = obj.CDF(obj.DiscreteX);
            obj.DiscretePDF = diff([0 obj.DiscreteCDF]);
            obj.DistType = 'd';  % return to discrete type
            TrimTables(obj,eps(0),1);
            SetBinEdges(obj);
        end

        function parmvals = ParmValues(obj,varargin)
            parmvals = [obj.OrderK];
            for iDist=1:obj.SampleSize
                parmvals = [parmvals obj.BasisRV{iDist}.ParmValues];
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % Improve here and next procedure by constraining Order <= SampleSize.
            Reals = NumTrans.GT2Real(1,Parms(1));
            ParmsUsed = 1;
            for iDist=1:obj.SampleSize
                Reals = [Reals obj.BasisRV{iDist}.ParmsToReals(Parms(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms))];
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(1,Reals(1));
            ParmsUsed = 1;
            for iDist=1:obj.SampleSize
                Parms = [Parms obj.BasisRV{iDist}.RealsToParms(Reals(ParmsUsed+1:ParmsUsed+obj.BasisRV{iDist}.NDistParms))];
                ParmsUsed = ParmsUsed + obj.BasisRV{iDist}.NDistParms;
            end
        end
        
        % function thisval=LegalValue(obj,X)
        %     thisval = LegalValue(obj.BasisRV{1},X);
        %     %if Not DLCreated Then CreateDL;
        %     %LegalValue = DiscreteList.LegalValue(X);
        % end
        
        % function thisval=NearestLegal(obj,X)
        %     thisval = NearestLegal(obj.BasisRV{1},X);
        %     %if Not DLCreated Then CreateDL;
        %     %NearestLegal = DiscreteList.NearestLegal(X);
        % end
        
        % function thisval=nIthValue(obj,I)
        %     thisval = nIthValue(obj.BasisRV{1},I);
        %     %if Not DLCreated Then CreateDL;
        %     %IthValue = DiscreteList.IthValue(I);
        % end
        
        function thiscdf=CDF(obj,X)
            if obj.DistType=='d'
                thiscdf = CDF@dDiscrete(obj,X);
                return;
            end
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            % Fx = zeros(ord.SampleSize+1,1);
            % One_Fx = zeros(ord.SampleSize+1,1);
            Fx = zeros(obj.SampleSize,1);
            for iel=1:numel(X)
                if InBounds(iel)
                    for I = 1:obj.SampleSize
                        Fx(I) = obj.BasisRV{I}.CDF(X(iel));
                    end
                    One_Fx = 1 - Fx;
                    % Now look at all possible combinations of Order BasisRV's
                    % less than X and the rest greater than X.  Compute
                    % the probability of each such combination and sum across those.
                    FxSum = 0;
                    for I = obj.CDFLoop1:obj.CDFLoop2
                        % Determine all possible subsets for use in computing CDFs.
                        % Improvement needed: It would be faster to do this in constructor for all I in loop and save.
                        SubsetsSmaller = nchoosek(1:obj.SampleSize,I);
                        NSubsetsSmaller = size(SubsetsSmaller,1);
                        for iSubset=1:NSubsetsSmaller
                            Indicators = zeros(obj.SampleSize,1);
                            Indicators(SubsetsSmaller(iSubset,:)) = 1;
                            One_FxTerm = 1;
                            FxTerm = 1;
                            for J = 1:obj.SampleSize
                                if Indicators(J)
                                    FxTerm = FxTerm * Fx(J);
                                else
                                    One_FxTerm = One_FxTerm * One_Fx(J);
                                end
                            end
                            FxSum = FxTerm * One_FxTerm + FxSum;
                        end
                    end
                    if obj.BigOrder
                        thiscdf(iel) = FxSum;
                    else
                        thiscdf(iel) = 1 - FxSum;
                    end
                end
            end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval=zeros(varargin{:});
            X = zeros(obj.SampleSize,1);
            % I really wonder how to vectorize this!
            for i=1:numel(thisval)
                for iSamp = 1:obj.SampleSize
                    X(iSamp) = Random(obj.BasisRV{iSamp});
                end
                X = sort(X);
                thisval(i) = X(obj.OrderK);
            end
        end
        
    end  % methods
    
end  % class Order

