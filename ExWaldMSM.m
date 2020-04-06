classdef ExWaldMSM < ExWald  % Old version superceded by simple convolution which seems much faster
    % ExWaldMSM(WaldMean, WaldSD, ExpMean): 3-parameter Ex-Wald distribution, the sum of a 2-parameter Wald, specified in terms of its mean & sd,
    % plus an exponential, specified in terms of its mean.  The Wald's sigma is 1.
    
    properties(SetAccess = protected)
        WaldMean, WaldSD, ExpMean
    end

    properties(SetAccess = public)
        MinWaldMean, MinWaldSD, MinExpMean
    end
    
    methods
        
        function obj=ExWaldMSM(varargin)
            obj=obj@ExWald;
            obj.FamilyName = 'ExWaldMSM';
            obj.NDistParms = 3;
            obj.ParmNames{1} = 'WaldMean';
            obj.ParmNames{2} = 'WaldSD';
            obj.ParmNames{3} = 'ExpMean';
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.MinWaldMean = 0.001; % These lower limits are to prevent numerical errors
            obj.MinWaldSD   = 0.001; % in parameter searching. .01 worked for all.
            obj.MinExpMean  = 0.001;
            obj.WaldBasis.InverseCDFTol = 0.0005;   % Tolerance on how close the CDF value must be to the target.
            obj.WaldBasish.InverseCDFTol = 0.0005;  % Tolerance on how close the CDF value must be to the target.
            % This would improve speed but I don't want to hard-code it in; user may call it after making the dist.
            % obj.SearchOptions.TolFun = 1e-4;
            % obj.SearchOptions.TolX = 1e-4;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExWaldMSM:Constructor', ...
                        'ExWaldMSM constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.WaldMean = newparmvalues(1);
            obj.WaldSD = newparmvalues(2);
            obj.ExpMean = newparmvalues(3);
            if obj.WaldMean<=obj.MinWaldMean
                error(['ExWaldMSM WaldMean must be > ' num2str(obj.MinWaldMean) '.']);
            end
            if obj.WaldSD<obj.MinWaldSD
                error(['ExWaldMSM WaldSD must be > ' num2str(obj.MinWaldSD) '.']);
            end
            if obj.ExpMean<obj.MinExpMean
                error(['ExWaldMSM ExpMean must be > ' num2str(obj.MinExpMean) '.']);
            end
            PassA = sqrt(obj.WaldMean)*obj.WaldMean/obj.WaldSD;
            PassMu = PassA / obj.WaldMean;
            PassRate = 1/obj.ExpMean;
            HoldNameBuilding = obj.NameBuilding;
            obj.NameBuilding = false;
            obj.NDistParms = 4;  % Temporarily for ExWald processing.
            ResetParms@ExWald(obj,[PassMu 1 PassA PassRate]);
            obj.NDistParms = 3;  % Return to the true value for ExWaldMSM.
            obj.NameBuilding = HoldNameBuilding;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newWaldMean = ifelse(ParmCodes(1)=='f',obj.WaldMean,1.05*obj.WaldMean);
            newWaldSD   = ifelse(ParmCodes(2)=='f',obj.WaldSD,  1.05*obj.WaldSD);
            newExpMean  = ifelse(ParmCodes(3)=='f',obj.ExpMean, 1.05*obj.ExpMean);
            obj.ResetParms([newWaldMean newWaldSD newExpMean]);
        end
        
        function []=ReInit(obj)
            HoldNameBuilding = obj.NameBuilding;
            obj.NameBuilding = false;
            ReInit@ExWald(obj);
            obj.NameBuilding = HoldNameBuilding;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(obj.MinWaldMean,Parms(1)) NumTrans.GT2Real(obj.MinWaldSD,Parms(2)) NumTrans.GT2Real(obj.MinExpMean,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(obj.MinWaldMean,Reals(1)) NumTrans.Real2GT(obj.MinWaldSD,Reals(2)) NumTrans.Real2GT(obj.MinExpMean,Reals(3))];
        end
        
    end  % methods
    
end  % class ExWaldMSM



