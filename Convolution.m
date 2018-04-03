classdef Convolution < dContinuous  % dEither
    % Convolution(BasisRV1,BasisRV2) creates a random variable that is
    %  the sum of two independent basis random variables,
    %  BasisRV1+BasisRV2
    
    % PDF only works for continuous distributions.
    % Should handle Convolutions of discrete RVs via ListRV mechanism.
    % When testing, remember that the distribution of the mean of 2 Uniform(0,1) is Triangular(0,1)
    
    properties(SetAccess = protected)
        BasisRV1, BasisRV2
    end
    
    methods
        
        function obj=Convolution(varargin)
            obj=obj@dContinuous('Convolution'); % dEither('Convolution');
            switch nargin
                case 0
                case 2
                    obj.BasisRV1 = varargin{1};
                    obj.BasisRV2 = varargin{2};
                    if (obj.BasisRV1.DistType=='c') && (obj.BasisRV2.DistType=='c')
                        obj.DistType = 'c';
                    else
                        assert(false,'Convolution can only handle continuous Basis distributions (so far)');
                    end
                    obj.NDistParms = obj.BasisRV1.NDistParms + obj.BasisRV2.NDistParms;
                    obj.DefaultParmCodes = [obj.BasisRV1.DefaultParmCodes obj.BasisRV2.DefaultParmCodes];
                    ResetParms(obj,[obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
                otherwise
                    ME = MException('Convolution:Constructor', ...
                        'Convolution constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function BuildMyName(obj)
            assert(obj.Initialized,UninitializedError(obj));
            obj.StringName = ['Convolution(' obj.BasisRV1.StringName ',' obj.BasisRV2.StringName ')'];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    thispdf(iel)=integral(@(x) obj.BasisRV1.PDF(X(iel)-x).*obj.BasisRV2.PDF(x),obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound);
                end
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,ParmCodes)
            %           assert(obj.Initialized,UninitializedError(obj));
            Reals = [obj.BasisRV1.ParmsToReals(Parms(1:obj.BasisRV1.NDistParms),ParmCodes(1:obj.BasisRV1.NDistParms)) ...
                obj.BasisRV2.ParmsToReals(Parms(obj.BasisRV1.NDistParms+1:end),ParmCodes(obj.BasisRV1.NDistParms+1:end)) ];
        end
        
        function Parms = RealsToParms(obj,Reals,ParmCodes)
            %           assert(obj.Initialized,UninitializedError(obj));
            Parms = [obj.BasisRV1.RealsToParms(Reals(1:obj.BasisRV1.NDistParms),ParmCodes(1:obj.BasisRV1.NDistParms)) ...
                obj.BasisRV2.RealsToParms(Reals(obj.BasisRV1.NDistParms+1:end),ParmCodes(obj.BasisRV1.NDistParms+1:end)) ];
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.BasisRV1.ResetParms(newparmvalues(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.ResetParms(newparmvalues(obj.BasisRV1.NDistParms+1:end));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            obj.BasisRV1.PerturbParms(ParmCodes(1:obj.BasisRV1.NDistParms));
            obj.BasisRV2.PerturbParms(ParmCodes(obj.BasisRV1.NDistParms+1:end));
            obj.ResetParms([obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.LowerBound = obj.BasisRV1.LowerBound + obj.BasisRV2.LowerBound;
            obj.UpperBound = obj.BasisRV1.UpperBound + obj.BasisRV2.UpperBound;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function parmvals = ParmValues(obj)
            % Return values of all parameters
            parmvals = [obj.BasisRV1.ParmValues obj.BasisRV2.ParmValues];
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Mean + obj.BasisRV2.Mean;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Variance + obj.BasisRV2.Variance;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Random(varargin{:}) + obj.BasisRV2.Random(varargin{:});
        end
        
        % function thisval=LegalValue(obj,X)
        %     thisval = LegalValue(obj.BasisRV{1},X);
        %     %if Not DLCreated Then CreateDL;
        %     %LegalValue = DiscreteList.LegalValue(X);
        % end
        
        % function thisval=NearestLegal(obj,X)
        %     thisval = NearestLegal(obj.BasisRV{1},X);
        %     %if Not DLCreated Then CreateDL;
        %     %NearestLegal = DiscreteList.NearestLegal(X);
        % end
        
        % function thisval=nIthValue(obj,I)
        %     thisval = nIthValue(obj.BasisRV{1},I);
        %     %if Not DLCreated Then CreateDL;
        %     %IthValue = DiscreteList.IthValue(I);
        % end
        
    end  % methods
    
end  % class Convolution

