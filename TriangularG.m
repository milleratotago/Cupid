classdef TriangularG < dContinuous
    % TriangularG(min,peak,max):  peak anywhere between min and max.
    % This distribution has serious problems during estimation because parameters cross over.  Try TriangularGCWP instead.
    
    properties(SetAccess = protected)
        min, peak, max,
        lowerarea, upperarea, peakheight,
        center, heightfac
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            % It is tempting to enforce the restrictions Parms(1) < Parms(2) < Parms(3), like this:
            % Reals = [Parms(1) NumTrans.GT2Real(Parms(1),Parms(2)) NumTrans.GT2Real(Parms(2),Parms(3))];
            % Unfortunately, this creates problems if you want to estimate some parms while holding other fixed.
            % It may be better to use a different parameterization to make the parameters independent.
            Reals = Parms;
        end
        
        function Parms = RealsToParms(Reals,~)
            %            temp = NumTrans.Real2GT(Reals(1),Reals(2));
            %            Parms = [Reals(1) NumTrans.Real2GT(Reals(1),Reals(2)) NumTrans.Real2GT(temp,Reals(3))];
            Parms = Reals;
        end
        
    end
    
    methods
        
        function obj=TriangularG(varargin)
            obj=obj@dContinuous('TriangularG');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('TriangularG:Constructor', ...
                        'TriangularG constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.min = newparmvalues(1);
            obj.peak = newparmvalues(2);
            obj.max = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Make the bounds a bit wider & move the peak down a little:
            OldLower = obj.min;
            OldUpper = obj.max;
            OldPeak = obj.peak;
            BoundShift = 0.051 * (OldUpper - OldLower);
            OldFrac = (OldPeak - OldLower) / (OldUpper - OldLower);
            NewLower = ifelse(ParmCodes(1)=='f',obj.LowerBound,OldLower-BoundShift);
            NewUpper = ifelse(ParmCodes(3)=='f',obj.UpperBound,OldUpper+BoundShift);
            NewPeak  = ifelse(ParmCodes(2)=='f',obj.peak,NewLower + 0.95*OldFrac*(NewUpper-NewLower));
            % if (NewLower > NewPeak) || (NewPeak > NewUpper)
            %     aOldLower = OldLower
            %     aOldUpper = OldUpper
            %     aOldPeak = OldPeak
            %     aBoundShift = BoundShift
            %     aOldFrac = OldFrac
            %     disp(['Perturb set: ' num2str(NewLower) ' ' num2str(NewPeak) ' ' num2str(NewUpper)]);
            %     pause;
            % end
            obj.ResetParms([NewLower NewPeak NewUpper]);
        end
        
        function []=ReInit(obj)
            if ~(obj.min<obj.peak)
                error(['TriangularG must satisfy min<peak but these are ' num2str(obj.min) ' and ' num2str(obj.peak)]);
            end
            if ~(obj.peak<obj.max)
                error(['TriangularG must satisfy peak<max but these are ' num2str(obj.peak) ' and ' num2str(obj.max)]);
            end
            obj.center = (obj.min + obj.max) / 2;
            obj.heightfac = 4 / (obj.max-obj.min)^2;
            obj.LowerBound = obj.min;
            obj.UpperBound = obj.max;
            obj.peakheight = 2 / (obj.max - obj.min);
            obj.lowerarea = 0.5 * obj.peakheight * (obj.peak - obj.min);
            obj.upperarea = 1 - obj.lowerarea;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            LowerPart = (X>=obj.LowerBound) & (X<=obj.peak);
            thispdf(LowerPart) = (X(LowerPart) - obj.min) / (obj.peak - obj.min) * obj.peakheight;
            UpperPart = (X<=obj.UpperBound) & (X>obj.peak);
            thispdf(UpperPart) = (obj.max - X(UpperPart)) / (obj.max - obj.peak) * obj.peakheight;
            % for i=1:numel(X)
            %     if (X(i) <= obj.min) || (X(i) >= obj.max)
            %     elseif X(i) == obj.peak
            %         thispdf(i) = obj.peakheight;
            %     elseif X(i) < obj.peak
            %         thispdf(i) = (X(i) - obj.min) / (obj.peak - obj.min) * obj.peakheight;
            %     else
            %         % if X(i) > obj.peak
            %         thispdf(i) = (obj.max - X(i)) / (obj.max - obj.peak) * obj.peakheight;
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            LowerPart = (X>=obj.LowerBound) & (X<=obj.peak);
            thiscdf(LowerPart) = 0.5* (X(LowerPart) - obj.min).^2 / (obj.peak - obj.min) * obj.peakheight;
            UpperPart = (X<=obj.UpperBound) & (X>obj.peak);
            thiscdf(UpperPart) = 1 - 0.5 * (obj.max - X(UpperPart)).^2 / (obj.max - obj.peak) * obj.peakheight;
            thiscdf(X>obj.UpperBound) = 1;
            % for i=1:numel(X)
            %     if X(i) <= obj.min
            %         thiscdf(i) = 0;
            %     elseif X(i) >= obj.max
            %         thiscdf(i) = 1;
            %     elseif X(i) == obj.peak
            %         thiscdf(i) = obj.lowerarea;
            %     elseif X(i) < obj.peak
            %         thiscdf(i) = 0.5 * (X(i) - obj.min)^2 / (obj.peak - obj.min) * obj.peakheight;
            %     else
            %         % X(i) > obj.peak
            %         thiscdf(i) = 1 - 0.5 * (obj.max - X(i))^2 / (obj.max - obj.peak) * obj.peakheight;
            %     end
            % end
        end
        
    end  % methods
    
end  % class TriangularG

