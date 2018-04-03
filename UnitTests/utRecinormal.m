classdef utRecinormal < utContinuous;
    
    properties (ClassSetupParameter)
        parm1mu       = struct( 'p_05',.05 , 'p_005',.005 , 'p_0065',.0065);
        parm2sigma    = struct( 'p_01',.002 , 'p_001',.001 , 'p_002',.002);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utRecinormal(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1mu,parm2sigma)
            % Computations specific to the Recinormal distribution.
            testCase.Dist = Recinormal(parm1mu,parm2sigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            testCase.RawMomentRelTol(2) = 0.02;
            testCase.CenMomentAbsTol(2) = 0.02;
            testCase.MGFh = 1.0E-8; % MGF gets big very fast so this must be quite small.
%           testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utRecinormal


