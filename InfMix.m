classdef InfMix < dContinuous  % dEither
    % InfMix(ParentDist,ParmDist,MixParmNum)  % Optional MixParmNum if ParentDist only has one parm?
    
    %{
This random variable implements a mixture distribution in which a parameter of one
distribution (ParentDist) varies randomly according to a second distribution (ParmDist).
Suppose, for example, that we generated random numbers, x_i, with the
following 3 steps for each random number:

 1  Choose a random number, r, uniformly distributed from 10 to 20.
 2  Create a normal distribution with mean 0 and standard deviation r.
 3  Choose a random number, x, from the normal distribution just constructed.

The distribution of all the x_i's is an infinite mixture distribution--that is,
a mixture of an infinite number of normal distributions, varying in their SDs.
In CUPID, this distribution would be called:
    InfMix(Normal(0,1),Uniform(10,20),2)
which is read roughly as "an infinite mixture distribution formed from
a normal distribution with mean 0 and a standard deviation 1 that
varies uniformly from 10 to 20."
Note that the value of 1 for the standard deviation of the normal is actually
meaningless, because the values of the standard deviation are really determined
by the Uniform distribution.
    %}
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        ParentDist   % The basic distribution whose parameter is being mixed
        ParmDist     % The distribution of the varying parameter
        MixParmNum   % The number of the parameter that is being mixed (e.g., if mixing the SD of a normal,
                     % this value would be 2 because the SD is the second parameter of the normal distribution).
        ParentParms  % Parameter values of the parent distribution.  The value of the MixParmNum'th parameter changes.
        TypicalParm  % An initial value for the parameter that is being mixed
        NBParent, NBParm
    end
    
    properties(SetAccess = public)  % These properties can be set by anyone.
        WeightedPDFAbsTol, WeightedPDFRelTol
        WeightedCDFAbsTol, WeightedCDFRelTol
        WeightedRawMomentAbsTol, WeightedRawMomentRelTol
    end
    
    methods
        
        function obj=InfMix(varargin)   % Constructor
            obj=obj@dContinuous('InfMix');  % Inherited constructor
            obj.WeightedPDFAbsTol = obj.IntegralPDFAbsTol;
            obj.WeightedPDFRelTol = obj.IntegralPDFRelTol;
            obj.WeightedCDFAbsTol = obj.IntegralCDFAbsTol;
            obj.WeightedCDFRelTol = obj.IntegralCDFAbsTol;
            obj.WeightedRawMomentAbsTol = obj.IntegralPDFXNAbsTol;
            obj.WeightedRawMomentRelTol = obj.IntegralPDFXNRelTol;
            obj.NBParent = zeros(10,1);
            obj.NBParm = zeros(10,1);
            switch nargin
                case 0
                case 3
                    Setup(obj,varargin(:));
                    ResetParms(obj,[obj.ParmValues]);
                otherwise
                    ME = MException('InfMix:Constructor', ...
                        'Illegal number of arguments passed to InfMix constructor.');
                    throw(ME);
            end
        end
        
        function Setup(obj,s)
            obj.ParentDist = s{1};
            obj.ParmDist = s{2};
            obj.MixParmNum = s{3};
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParms(obj);
            obj.ParentDist.ResetParms(newparmvalues(1:obj.ParentDist.NDistParms));
            obj.ParmDist.ResetParms(newparmvalues(obj.ParentDist.NDistParms+1:end));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            obj.ParentDist.PerturbParms(ParmCodes(1:obj.ParentDist.NDistParms));
            obj.ParmDist.PerturbParms(ParmCodes(obj.ParentDist.NDistParms+1:end));
            obj.ResetParms([obj.ParentDist.ParmValues obj.ParmDist.ParmValues]);
        end
        
        function []=ReInit(obj)
            % Re-initialize after parameters have been reset.
            assert((obj.MixParmNum>=1)&(obj.MixParmNum<=obj.ParentDist.NDistParms)&iswholenumber(obj.MixParmNum), ...
                'InfMix.MixParmNum must be an integer between 1 and the number of parameters in the parent distribution.');
            obj.ParmTypes = [obj.ParentDist.DefaultParmCodes obj.ParmDist.DefaultParmCodes];
            obj.DefaultParmCodes = [obj.ParentDist.DefaultParmCodes obj.ParmDist.DefaultParmCodes];
            obj.NDistParms = obj.ParentDist.NDistParms + obj.ParmDist.NDistParms;
            % obj.DistType = ParentDist.DistType;
            
            % The following lines implement a rather crude search to find the lower and
            % upper bounds of this mixture distribution.  Essentially, this algorithm
            % assumes that the upper & lower bounds of the mixture distribution occur when
            % ParmDist is at its upper & lower bounds, not necessarily in that order.
            % Note that it would be easy to define this as a separate function and let
            % the user over-ride it to define the boundss via a user function, but I
            % cannot think of an example where anyone would want to do that.
            obj.ParentParms = obj.ParentDist.ParmValues;
            obj.TypicalParm = obj.ParentParms(obj.MixParmNum);
            obj.ParentParms(obj.MixParmNum) = obj.ParmDist.LowerBound;
            obj.ParentDist.ResetParms(obj.ParentParms);
            Low1 = obj.ParentDist.LowerBound;
            Upp1 = obj.ParentDist.UpperBound;
            obj.ParentParms(obj.MixParmNum) = obj.ParmDist.UpperBound;
            obj.ParentDist.ResetParms(obj.ParentParms);
            Low2 = obj.ParentDist.LowerBound;
            Upp2 = obj.ParentDist.UpperBound;
            if Low1 < Low2
                obj.LowerBound = Low1;
            else
                obj.LowerBound = Low2;
            end
            if Upp1 > Upp2
                obj.UpperBound = Upp1;
            else
                obj.UpperBound = Upp2;
            end
            
            obj.ParentParms(obj.MixParmNum) = obj.TypicalParm;
            obj.ParentDist.ResetParms(obj.ParentParms);
            % dEither: This only works if the ParentDist distribution has all of
            % its values with the typical value of ParmDist.
            % if DistType = Discrete Then NValues = ParentDist.NValues;
            
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function PushAndStopNameBuilding(obj)
            PushAndStopNameBuilding@dGeneric(obj);
            obj.ParentDist.NameBuilding = false;
            obj.ParmDist.NameBuilding = false;
        end

        function PopNameBuilding(obj)
            obj.ParentDist.NameBuilding = obj.NBParent(obj.NBNStored);
            obj.ParmDist.NameBuilding = obj.NBParm(obj.NBNStored);
            PopNameBuilding@dGeneric(obj);  % Decrements NBStored so this must go last
        end

        function BuildMyName(obj)
            obj.StringName = ['InfMix(' obj.ParentDist.StringName ',' obj.ParmDist.StringName ',' num2str(obj.MixParmNum) ')'];
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            % Convert the current parameter values to a list of reals in the range (-inf,inf) for fminsearch to adjust.
            Reals = [obj.ParentDist.ParmsToReals(Parms(1:obj.ParentDist.NDistParms)) obj.ParmDist.ParmsToReals(Parms(obj.ParentDist.NDistParms+1:end))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            % Convert fminsearch's reals from the range (-inf,inf) into a list of legal parameter values.
            Parms = [obj.ParentDist.RealsToParms(Reals(1:obj.ParentDist.NDistParms)) obj.ParmDist.RealsToParms(Reals(obj.ParentDist.NDistParms+1:end))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.ParentDist.ParmValues obj.ParmDist.ParmValues];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if obj.DistType == 'd'
                error('Discrete InfMix not implemented.');
                % thispdf(InBounds) = zeros(size(X));
            else
                obj.PushAndStopNameBuilding;
                for iel=1:numel(X)
                    if InBounds(iel)
                        thispdf(iel) = integral(@(x) WeightedPDF(obj,X(iel),x),obj.ParmDist.LowerBound,obj.ParmDist.UpperBound,'AbsTol',obj.WeightedPDFAbsTol,'RelTol',obj.WeightedPDFRelTol);
                    end
                end
                obj.PopNameBuilding;
            end
        end
        
        function thisval=WeightedPDF(obj,X,parmval)
            thisval = zeros(size(parmval));
            for iel=1:numel(parmval)
                Pr = obj.ParmDist.PDF(parmval(iel));
                obj.ParentParms(obj.MixParmNum) = parmval(iel);
                obj.ParentDist.ResetParms(obj.ParentParms);
                thisval(iel) = obj.ParentDist.PDF(X) * Pr;
            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            if obj.DistType == 'd'
                error('Discrete InfMix not implemented.');
                % thiscdf(InBounds) = zeros(size(X));
            else
                obj.PushAndStopNameBuilding;
                for iel=1:numel(X)
                    if InBounds(iel)
                        thiscdf(iel) = integral(@(x) WeightedCDF(obj,X(iel),x),obj.ParmDist.LowerBound,obj.ParmDist.UpperBound,'AbsTol',obj.WeightedCDFAbsTol,'RelTol',obj.WeightedCDFRelTol);
                    end
                end
                obj.PopNameBuilding;
            end
        end
        
        function thisval=WeightedCDF(obj,X,parmval)
            thisval = zeros(size(parmval));
            for iel=1:numel(parmval)
                Pr = obj.ParmDist.PDF(parmval(iel));
                obj.ParentParms(obj.MixParmNum) = parmval(iel);
                obj.ParentDist.ResetParms(obj.ParentParms);
                thisval(iel) = obj.ParentDist.CDF(X) * Pr;
            end
        end
        
        function thisval=RawMoment(obj,I)
            if obj.DistType == 'd'
                error('Discrete InfMix not implemented.');
                % thisval = 0;
            else
                if I >= 1
                    AbsTol = obj.WeightedRawMomentAbsTol(I);
                    RelTol = obj.WeightedRawMomentRelTol(I);
                else
                    AbsTol = obj.WeightedPDFAbsTol;
                    RelTol = obj.WeightedPDFRelTol;
                end
                obj.PushAndStopNameBuilding;
                thisval = integral(@(x) WeightedMoment(obj,I,x),obj.ParmDist.LowerBound,obj.ParmDist.UpperBound,'AbsTol',AbsTol,'RelTol',RelTol);
                obj.PopNameBuilding;
            end
        end
        
        function thisval=WeightedMoment(obj,I,parmval)
            thisval = zeros(size(parmval));
            for iel=1:numel(parmval)
                Pr = obj.ParmDist.PDF(parmval(iel));
                obj.ParentParms(obj.MixParmNum) = parmval(iel);
                obj.ParentDist.ResetParms(obj.ParentParms);
                thisval(iel) = obj.ParentDist.RawMoment(I) * Pr;
            end
        end
        
        function thisval=CenMoment(obj,I)
            thisval=CenMomentFromRawMoment(obj,I);
        end

        function thisval=Random(obj,varargin)
            thisval=zeros(varargin{:});
            randparms = obj.ParmDist.Random(varargin{:});
            obj.PushAndStopNameBuilding;
            for iel=1:numel(thisval)
                obj.ParentParms(obj.MixParmNum) = randparms(iel);
                obj.ParentDist.ResetParms(obj.ParentParms);
                thisval(iel) = obj.ParentDist.Random;
            end
            obj.PopNameBuilding;
        end
        
    end  % methods
    
end  % object InfMix

%{

%*********InfMix Distributions:*******************************
The following object implements InfMix random variables.  They are derived
from a ParentDist distribution, but the MixParm'th parameter of the ParentDist
distribution is distributed according to the ParmDist distribution.
}
TInfMixRV = Class(GenericRV)
   IMMom : TInfMixMoment;
   IMPDF : TInfMixPDF;
   IMCDF : TInfMixCDF;
   TypicalParm;
function thisval=Init(PassMixParm );
    % The ParentDist and ParmDist distributions must ALREADY have been initialized
      before Init is called!
function thisval=LegalValue(X);
function thisval=NearestLegal(X);
function thisval=IthValue(I );   % Returns Ith discrete value in sequence.


%*********InfMix Random Variables:****************************}


function thisval=TInfMixRV.LegalValue(X);
LegalValue = ParentDist.LegalValue(X);
end

function thisval=TInfMixRV.NearestLegal(X);
NearestLegal = ParentDist.NearestLegal(X);
end

function thisval=TInfMixRV.IthValue(I );
% NOTE: This only works if the ParentDist distribution has all of
% its values with the typical value of ParmDist.
ParentDist.ResetParm(MixParm,TypicalParm);
IthValue = ParentDist.IthValue(I);
end

function thisval=TInfMixRV.PDF(x);
Var This, OneVal, OnePr, OnePDF; HoldNB;
    I ;
HoldNB = NameBuilding;
NameBuilding = false;
This = 0;
Case ParmDist.DistType of
 Continuous : Begin
              IMPDF.XVal = X;
              if Not IMPDF.Integral(ParmDist.obj.LowerBound,ParmDist.obj.UpperBound,IntegralPrecisAbs,
                 IntegralPrecisRel,IntMinSteps,This) Then warning(
                 'TInfMixRV.PDF: Integral did not converge.');
              end
 Discrete : for I = 1:ParmDist.NValues Do Begin
               OneVal = ParmDist.IthValue(I);
               OnePr = ParmDist.PDF(OneVal);
               ParentDist.ResetParm(MixParm,OneVal);
               Onethispdf = ParentDist.PDF(X);
               This = This + OnePr*OnePDF;
               end
 else warning('Cannot compute InfMix PDF for mixed parameter distributions.');
 end
thispdf = This;
NameBuilding = HoldNB;
end

function thisval=TInfMixRV.CDF(x);
Var ThisC, ThisD, OneVal, OnePr, OneCDF; HoldNB;
    I ;
HoldNB = NameBuilding;
NameBuilding = false;
ThisC = 0;
ThisD = 0;
if ParmDist.DistType in (Continuous,Mixed)
   IMCDF.XVal = X;
   if Not IMCDF.Integral(ParmDist.obj.LowerBound,ParmDist.obj.UpperBound,IntegralPrecisAbs,
    IntegralPrecisRel,IntMinSteps,ThisC) Then warning(
     'TInfMixRV.CDF: Integral did not converge.');
   end
if ParmDist.DistType in (Discrete,Mixed)
   for I = 1:ParmDist.NValues Do Begin
      OneVal = ParmDist.IthValue(I);
      OnePr = ParmDist.PDF(OneVal);
      ParentDist.ResetParm(MixParm,OneVal);
      Onethiscdf = ParentDist.CDF(X);
      ThisD = ThisD + OnePr*OneCDF;
      end
   end
thiscdf = ThisC + ThisD;
NameBuilding = HoldNB;
end

function thisval=TInfMixRV.RawMoment(I );
Var ThisC, ThisD, OneVal, OnePr, OneMom; HoldNB;
    J ;
% Writeln('Entered TInfMixRV.RawMoment.');
HoldNB = NameBuilding;
NameBuilding = false;
ThisC = 0;
ThisD = 0;
if ParmDist.DistType in (Continuous,Mixed)
   IMMom.IMom = I;
   if Not IMMom.Integral(ParmDist.obj.LowerBound,ParmDist.obj.UpperBound,IntegralPrecisAbs,
    IntegralPrecisRel,IntMinSteps,ThisC) Then warning(
     'TInfMixRV.RawMoment: Integral did not converge.');
%  Writeln('TInfMixRV.RawMoment computed continuous Moment ',I,' as ',ThisC);
   end
if ParmDist.DistType in (Discrete,Mixed)
   for J = 1:ParmDist.NValues Do Begin
      OneVal = ParmDist.IthValue(J);
      OnePr = ParmDist.PDF(OneVal);
      ParentDist.ResetParm(MixParm,OneVal);
      OneMom = ParentDist.RawMoment(I);
      ThisD = ThisD + OnePr*OneMom;
%     Writeln('TInfMixRV.RawMoment computed discrete Moment ',I,' as ',ThisD);
      end
   end
RawMoment = ThisC + ThisD;
NameBuilding = HoldNB;
end

%}
