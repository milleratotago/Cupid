classdef Quantal < dContinuous
    % This is a distribution derived from quantal fluctuations in psychophysics.
    % The parameter is a positive integer.
    
    properties(SetAccess = protected)
        Threshold, TM1
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
        
        function obj=Quantal(varargin)
            obj=obj@dContinuous('Quantal');
            obj.NDistParms = 1;
            obj.ParmTypes = 'i';
            obj.DefaultParmCodes = 'i';
            switch nargin
                case 0
                case 1
                    ResetParms(obj,varargin{1});
                otherwise
                    ME = MException('Quantal:Constructor', ...
                        'Quantal constructor needs 0 or 1 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.Threshold = VerifyIntegerGE(obj,1,newparmvalues(1));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            if ~(ParmCodes(1)=='f')
                obj.ResetParms(obj.Threshold+1);
            end
        end
        
        function []=ReInit(obj)
            obj.TM1 = obj.Threshold - 1;
            obj.LowerBound = 0;
            obj.UpperBound = 10 * obj.Threshold;
            obj.Initialized = true;
            while CDF(obj,obj.UpperBound) < obj.CDFNearlyOne
                obj.UpperBound = obj.UpperBound + 5;
            end
            % obj.UpperBound = InverseCDF(obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        %        function thispdf=PDF(obj,X)
        %            thispdf=PDFfromCDF(obj,X);
        %        end
        %
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            PSum = zeros(size(X));
            XFac = 1;
            for I = 0:obj.TM1
                PSum(InBounds) = PSum(InBounds) + X(InBounds).^I / XFac .* exp(-X(InBounds));
                XFac = XFac * (I+1);
            end
            thiscdf(InBounds) = 1 - PSum((InBounds));
        end
        
    end  % methods
    
end  % class Quantal
