classdef ProdUni01 < dContinuous
    % ProdUni01(K) is the distribution of the product of K independent Uniform(0,1) random variables.
    % See https://math.stackexchange.com/questions/659254/product-distribution-of-two-uniform-distribution-what-about-3-or-more
    % To get the distribution of the product of Uniform(0,p) random variables, use MultTrans(ProdUni01(K),p^K), for p<1
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        K
        Km1, Km1fac
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = NumTrans.GT2Real(1,Parms(1));
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = NumTrans.Real2GT(1,Reals(1));
        end
        
    end
    
    methods
        
        function obj=ProdUni01(varargin)
            obj=obj@dContinuous('ProdUni01');
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            obj.NDistParms = 1;
            obj.LowerBound = eps;
            obj.UpperBound = 1 - eps;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ProdUni01:Constructor', ...
                        'ProdUni01 constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.K = round(newparmvalues(1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newK = ifelse(ParmCodes(1)=='f', obj.K,obj.K+1);
            obj.ResetParms(newK);
        end
        
        function []=ReInit(obj)
            if (obj.K<1) || (floor(obj.K) ~= obj.K)
                error('ProdUni01 K must be an integer > 0.');
            end
            obj.Km1 = obj.K - 1;
            obj.Km1fac = factorial(obj.Km1);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = (-1)^obj.Km1 * log(X(InBounds)).^obj.Km1 / obj.Km1fac;
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = rand(varargin{:},obj.K);
            thisval = prod(thisval,ndims(thisval));
        end
        
    end  % methods
    
end  % class ProdUni01
