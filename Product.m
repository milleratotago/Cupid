classdef Product < dContinuous  % dEither
    % Product(BasisRV1,BasisRV2) creates a random variable that is
    %  the Product of the two independent POSITIVE basis random variables,
    %  BasisRV1*BasisRV2
    
    properties(SetAccess = protected)
        BasisRV1, BasisRV2
        CDFRelTol, CDFAbsTol
    end
    
    methods
        
        function obj=Product(varargin)
            obj=obj@dContinuous('Product');
            obj.CDFRelTol = 1e-3;
            obj.CDFAbsTol = 1e-4;
            switch nargin
                case 0
                case 2
                    obj.BasisRV1 = varargin{1};
                    obj.BasisRV2 = varargin{2};
                    if (obj.BasisRV1.DistType=='c') && (obj.BasisRV2.DistType=='c')
                        obj.DistType = 'c';
                    else
                        assert(false,'Product can only handle continuous Basis distributions (so far)');
                    end
                    obj.NDistParms = obj.BasisRV1.NDistParms + obj.BasisRV2.NDistParms;
                    obj.DefaultParmCodes = [obj.BasisRV1.DefaultParmCodes obj.BasisRV2.DefaultParmCodes];
                    ResetParms(obj,[obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
                otherwise
                    ME = MException('Product:Constructor', ...
                        'Product constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = ['Product(' obj.BasisRV1.StringName ',' obj.BasisRV2.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.BasisRV1.ResetParms(newparmvalues(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.ResetParms(newparmvalues(obj.BasisRV1.NDistParms+1:end));
            assert((obj.BasisRV1.LowerBound>0)&&(obj.BasisRV2.LowerBound>0),'Product limited to positive RVs.');
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV1.PerturbParms(ParmCodes(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.PerturbParms(ParmCodes(obj.BasisRV1.NDistParms+1:end));
            obj.ResetParms([obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = obj.BasisRV1.LowerBound * obj.BasisRV2.LowerBound;
            obj.UpperBound = obj.BasisRV1.UpperBound * obj.BasisRV2.UpperBound;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    fun = @(x) obj.BasisRV1.CDF(X(iel)./x).*obj.BasisRV2.PDF(x);
                    thiscdf(iel)=integral(fun,obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,'RelTol', obj.CDFRelTol, 'AbsTol', obj.CDFAbsTol);
                end
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
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Random(varargin{:}) .* obj.BasisRV2.Random(varargin{:});
        end
        
    end  % methods
    
end  % class Product


