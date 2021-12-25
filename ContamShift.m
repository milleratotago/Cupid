classdef ContamShift < Mixture
    % ContamShift(SingleBasisRV,pContam,ShiftRV) creates a random variable that is
    % a mixture of the indicated basis RV and a shifted version of that RT,
    % with the size of the shift coming from the ShiftRV distribution.
    % pContam is the probability that the shift is added (i.e., that the SingleBasisRV is "contaminated").
    
    properties(SetAccess = protected)
        SingleBasisRV % The original uncontaminated RV.
        ShiftRV       % The distribution of increments added to the SingleBasisRV by contamination.
        ContamRV      % The convolution distribution of SingleBasisRV+ShiftRV
        pContam       % The probability of contamination.
        SBparms       % Positions within parmlists corresponding to SingleBasisRV
        Pparm         % Position within parmlists corresponding to pContam
        Shiftparms    % Positions within parmlists corresponding to ShiftRV
    end
    
    methods
        
        function obj=ContamShift(varargin)
            obj=obj@Mixture;
            obj.FamilyName = 'ContamShift';
            obj.TrimBounds = true;  % to avoid warnings about impossible values where PDF is 0 for observations in bounds.
            switch nargin
                case 3
                    Setup(obj,varargin(:));
                    ResetParms(obj,[obj.ParmValues]);
                otherwise
                    ME = MException('ContamShift:Constructor', ...
                        'ContamShift constructor needs 0 or 3+ arguments.');
                    throw(ME);
            end
        end

        function Setup(obj,s)
            % Just take the 3 parameters of ContamShift and construct the parameters
            % needed for Mixture's Setup.
            obj.SingleBasisRV = s{1};
            obj.pContam = s{2};
            obj.ShiftRV = s{3};
            obj.SBparms = 1:obj.SingleBasisRV.NDistParms;
            obj.Pparm = obj.SingleBasisRV.NDistParms + 1;
            obj.Shiftparms = (obj.Pparm+1):(obj.SingleBasisRV.NDistParms+1+obj.ShiftRV.NDistParms);
            obj.ContamRV = Convolution(obj.SingleBasisRV,obj.ShiftRV);
            obj.DefaultParmCodes = [obj.SingleBasisRV.DefaultParmCodes 'r' obj.ShiftRV.DefaultParmCodes];
            obj.NDistParms = numel(obj.DefaultParmCodes);
            sCell = {1-obj.pContam, obj.SingleBasisRV, obj.ContamRV};
            Setup@Mixture(obj,sCell);
        end % Setup
        
        function BuildMyName(obj)
            obj.StringName = ['ContamShift(' obj.SingleBasisRV.StringName ',' num2str(obj.pContam) ',' obj.ShiftRV.StringName ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            % The newparmvalues are in the order SingleBasisRV parms, pContam, shiftRV parms.
            % Mixture wants 1-pContam, SingleBasisRV parms, SinglBasisRV parms again for convolu, shift parms
            SingleBasisParms = newparmvalues(obj.SBparms);
            obj.pContam = newparmvalues(obj.Pparm);
            ShiftParms = newparmvalues(obj.Shiftparms);
            mixtureNewparmvalues = [1-obj.pContam, SingleBasisParms, SingleBasisParms, ShiftParms];
            ResetParms@Mixture(obj,mixtureNewparmvalues);
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
            obj.ShiftRV.PerturbParms(ParmCodes(obj.Shiftparms));
            ResetParms(obj,[obj.ParmValues]);
 %            obj.ReInit;
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [obj.SingleBasisRV.ParmsToReals(Parms(obj.SBparms)) NumTrans.Bounded2Real(0,1,Parms(obj.Pparm)) obj.ShiftRV.ParmsToReals(Parms(obj.Shiftparms))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [obj.SingleBasisRV.RealsToParms(Reals(obj.SBparms)) NumTrans.Real2Bounded(0,1,Reals(obj.Pparm)) obj.ShiftRV.RealsToParms(Reals(obj.Shiftparms))];
        end
        
        function parmvals = ParmValues(obj)
            parmvals = [obj.SingleBasisRV.ParmValues, obj.pContam, obj.ShiftRV.ParmValues];
        end
        
    end  % methods
    
end  % class ContamShift




