classdef Rosin < Weibull
    % This is the Weibull distribution with origin=0
    % Rosin-Rammler particle size distribution (Dm, Power), arises in analyses of particle sizes.
    % https://en.wikipedia.org/wiki/Particle-size_distribution
    
    methods
        
        function obj=Rosin(varargin)
            obj=obj@Weibull;
            obj.FamilyName = 'Rosin';
            obj.ParmTypes = 'rrf';
            obj.DefaultParmCodes = 'rrf';
            obj.origin = 0;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Rosin:Constructor', ...
                        'Rosin constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            obj.StringName = [obj.FamilyName '(' num2str(obj.scale) ',' num2str(obj.power) ')'];
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);  % NWJEFF: Was ,[newparmvalues 0]
            obj.scale = newparmvalues(1);
            obj.power = newparmvalues(2);
            ReInit(obj);
        end
        
    end  % methods
    
end  % class Rosin

