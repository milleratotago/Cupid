classdef TriangularGCWP < TriangularG
    % TriangularGCWP(center,width,peakprop)
    
    properties(SetAccess = protected)
        width, peakprop
    end

    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [ Parms(1) NumTrans.GT2Real(0,Parms(2)) NumTrans.Bounded2Real(0,1,Parms(3)) ];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [ Reals(1) NumTrans.Real2GT(0,Reals(2)) NumTrans.Real2Bounded(0,1,Reals(3)) ];
        end
        
    end
    
    methods
        
        function obj=TriangularGCWP(varargin)
            obj=obj@TriangularG;
            obj.FamilyName = 'TriangularGCWP';
            obj.ParmNames{1} = 'center';
            obj.ParmNames{2} = 'width';
            obj.ParmNames{3} = 'peakprop';
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('TriangularGCWP:Constructor', ...
                        'TriangularGCWP constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.center = newparmvalues(1);
            obj.width = newparmvalues(2);
            obj.peakprop = newparmvalues(3);
            assert((obj.peakprop>0)&&(obj.peakprop<1),'TriangularGCWP Peakprop parameter must be between 0 and 1.');
            obj.min = obj.center - obj.width/2;
            obj.max = obj.center + obj.width/2;
            obj.peak = obj.min + (obj.max - obj.min) * obj.peakprop;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Make the bounds a bit wider & move the peakprop down a little:
            Newcenter = ifelse(ParmCodes(1)=='f',obj.center,obj.center+0.1*obj.width);
            Newwidth = ifelse(ParmCodes(2)=='f',obj.width,1.1*obj.width);
            NewPeakprop  = ifelse(ParmCodes(3)=='f',obj.peakprop,0.95*obj.peakprop);
            obj.ResetParms([Newcenter Newwidth NewPeakprop]);
        end
        
    end  % methods
    
end  % class TriangularGCWP

