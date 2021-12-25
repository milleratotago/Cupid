classdef ContamStretch < Mixture  % NEWJEFF: Lots of code duplication with ContamShift
    % ContamStretch(SingleBasisRV,pContam,StretchRV) creates a random variable that is
    % a mixture of the indicated basis RV and a stretched version of that RT,
    % with stretching by a multiplicative factor from the StretchRV distribution
    % (in practice I expect StretchRV.LowerBound >= 1).
    % pContam is the probability of stretching (i.e., that the SingleBasisRV is "contaminated").
    
    properties(SetAccess = protected)
        SingleBasisRV % The original uncontaminated RV.
        StretchRV     % The distribution of stretch factors by which SingleBasisRV is multiplied when there is contamination.
        ContamRV      % The product distribution of SingleBasisRV*StretchRV
        pContam       % The probability of contamination.
        SBparms       % Positions within parmlists corresponding to SingleBasisRV
        Pparm         % Position within parmlists corresponding to pContam
        Stretchparms  % Positions within parmlists corresponding to StretchRV
    end
    
    methods
        
        function obj=ContamStretch(varargin)
            obj=obj@Mixture;
            obj.FamilyName = 'ContamStretch';
            switch nargin
                case 3
                    Setup(obj,varargin(:));
                    ResetParms(obj,[obj.ParmValues]);
                otherwise
                    ME = MException('ContamStretch:Constructor', ...
                        'ContamStretch constructor needs 0 or 3+ arguments.');
                    throw(ME);
            end
        end

        function Setup(obj,s)
            % Just take the 3 parameters of ContamStretch and construct the parameters
            % needed for Mixture's Setup.
            obj.SingleBasisRV = s{1};
            obj.pContam = s{2};
            obj.StretchRV = s{3};
            obj.SBparms = 1:obj.SingleBasisRV.NDistParms;
            obj.Pparm = obj.SingleBasisRV.NDistParms + 1;
            obj.Stretchparms = (obj.Pparm+1):(obj.SingleBasisRV.NDistParms+1+obj.StretchRV.NDistParms);
            obj.ContamRV = Product(obj.SingleBasisRV,obj.StretchRV);
            obj.DefaultParmCodes = [obj.SingleBasisRV.DefaultParmCodes 'r' obj.StretchRV.DefaultParmCodes];
            obj.NDistParms = numel(obj.DefaultParmCodes);
            sCell = {1-obj.pContam, obj.SingleBasisRV, obj.ContamRV};
            Setup@Mixture(obj,sCell);
            obj.CDFNearlyOne = 0.99999;
        end % Setup
        
        function BuildMyName(obj)
            obj.StringName = ['ContamStretch(' obj.SingleBasisRV.StringName ',' num2str(obj.pContam) ',' obj.StretchRV.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            % The newparmvalues are in the order SingleBasisRV parms, pContam, StretchRV parms.
            % Mixture wants 1-pContam, SingleBasisRV parms, SinglBasisRV parms again for convolu, shift parms
            SingleBasisParms = newparmvalues(obj.SBparms);
            obj.pContam = newparmvalues(obj.Pparm);
            StretchParms = newparmvalues(obj.Stretchparms);
            mixtureNewparmvalues = [1-obj.pContam, SingleBasisParms, SingleBasisParms, StretchParms];
            ResetParms@Mixture(obj,mixtureNewparmvalues);
            obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Make sure that the perturbation of the SingleBasisRV parms is carried over to the convolution.
            obj.SingleBasisRV.PerturbParms(ParmCodes(obj.SBparms));
            if ParmCodes(obj.Pparm) ~= 'f'
                if obj.pContam > 0.5
                    obj.pContam = obj.pContam - rand/100;
                else
                    obj.pContam = obj.pContam + rand/100;
                end
            end
            obj.StretchRV.PerturbParms(ParmCodes(obj.Stretchparms));
            ResetParms(obj,[obj.ParmValues]);
 %            obj.ReInit;
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.SingleBasisRV.ParmsToReals(Parms(obj.SBparms)) NumTrans.Bounded2Real(0,1,Parms(obj.Pparm)) obj.StretchRV.ParmsToReals(Parms(obj.Stretchparms))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.SingleBasisRV.RealsToParms(Reals(obj.SBparms)) NumTrans.Real2Bounded(0,1,Reals(obj.Pparm)) obj.StretchRV.RealsToParms(Reals(obj.Stretchparms))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.SingleBasisRV.ParmValues, obj.pContam, obj.StretchRV.ParmValues];
        end
        
    end  % methods
    
end  % class ContamStretch




