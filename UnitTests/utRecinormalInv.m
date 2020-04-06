classdef utRecinormalInv < utContinuous;
    
    properties (ClassSetupParameter)
        parm1mu       = struct( 'p_05' ,1/0.05  , 'p_005',1/0.005 , 'p_0065' ,1/0.0065);
        parm2sigma    = struct( 'p_002',1/0.002 , 'p_001',1/0.001 , 'p_00022',1/0.00022);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utRecinormalInv(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1mu,parm2sigma)
            % Computations specific to the RecinormalInv distribution.
            testCase.Dist = RecinormalInv(parm1mu,parm2sigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.RawMomentRelTol(2) = 0.02;
            testCase.CenMomentAbsTol(2) = 0.02;
            testCase.MGFh = 1.0E-8; % MGF gets big very fast so this must be quite small.
%           testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            testCase.SkipAllEst = true;  % Too slow because it has to estimates underlying mu and sigma for each finalmu, finalsigma

            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utRecinormalInv


