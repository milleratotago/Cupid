classdef Weibull2 < Weibull
    % Weibull2(scale>0,power>0)  origin=0
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.Bounded2Real(Weibull.minPower,Weibull.maxPower,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2Bounded(Weibull.minPower,Weibull.maxPower,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=Weibull2(varargin)
            obj=obj@Weibull;
            obj.origin = 0;
            obj.FamilyName = 'Weibull2';
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Weibull2:Constructor', ...
                        'Weibull2 constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.scale = newparmvalues(1);
            obj.power = newparmvalues(2);
            ReInit(obj);
        end
        
       function parms = ParmValues(obj)
            parms = [obj.scale, obj.power];
        end
        
       function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newscale  = ifelse(ParmCodes(1)=='f', obj.scale, 1.05*obj.scale);
            newpower  = ifelse(ParmCodes(2)=='f', obj.power, 0.95*obj.power);
            obj.ResetParms([newscale newpower]);
        end
        
        function []=ReInit(obj)
            assert(obj.scale>0,'Weibull2 scale must be > 0.');
            assert(obj.power>0,'Weibull2 power must be > 0.');
            obj.invpower = 1 / obj.power;
            obj.Initialized = true;
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parms = StartParmsMLE(obj,X)
            obsmedian = median(X);
            obs90pct = prctile(X,90);
            syms shape tttt
            estscale = obsmedian;
            F(tttt) = 1 - exp( -(tttt/estscale)^shape);
            vpsoln = vpasolve(F(obs90pct)==0.90,shape);
            estshape = double(vpsoln);
            parms = [estscale, estshape];
        end
        
    end  % methods
    
end  % class Weibull2

