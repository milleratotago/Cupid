classdef Frechet2 < Frechet
    % Frechet2 special case of Frechet distribution with minval = 0
    % Same as the ExtrVal2L distribution.
    
    methods(Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)), NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)), NumTrans.Real2GT(eps,Reals(2))];
        end
        
    end % methods(Static)
    
    methods
        
        function obj=Frechet2(varargin)
            obj=obj@Frechet;
            obj.FamilyName = 'Frechet2';
            obj.NDistParms = 2;
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.CDFNearlyOne = 0.999999;
            obj.minval = 0;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{1},varargin{2}]);
                otherwise
                    ME = MException('Frechet2:Constructor', ...
                        'Frechet2 constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.shape = newparmvalues(1);
            obj.scale = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newShape = ifelse(ParmCodes(1)=='f', obj.shape, 0.9*obj.shape);
            newScale = ifelse(ParmCodes(2)=='f', obj.scale, 0.9*obj.scale);
            obj.ResetParms([newShape, newScale]);
        end
        
        function []=ReInit(obj)
            assert(obj.shape>0,'Frechet2 shape parameter must be > 0.');
            assert(obj.scale>0,'Frechet2 scale parameter must be > 0.');
            obj.LowerBound = obj.minval+eps;
            obj.UpperBound = Frechet.frechicdf(obj.CDFNearlyOne,obj.shape,obj.scale,obj.minval);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parms = StartParmsMLE(obj,X)
            % Based on the following quartile functions from Wikipedia:
            % assuming m == 0
            % q1 = s * log(4)^(-1/alpha)
            % q3 = s * log(4/3)^(-1/alpha)
            obs = double(prctile(X,[25 75]));
            % Here is the array of functions, all of which are to be zero'ed
            F = @(p) [...
                ( obs(1)/p(2) )^(-p(1)) - log(4); ...
                ( obs(2)/p(2) )^(-p(1)) - log(4/3) ...
                ];
            starts = mean(X);
            x0 = double([3 starts]);
            parms = fsolve(F,x0,obj.fsolveoptions);
        end
        
    end  % methods
    
end  % class Frechet2

