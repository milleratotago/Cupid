classdef utRayleigh < utContinuous
    
    properties (ClassSetupParameter)
        parmSigma  = struct( 'p_1',.1 , 'p1',1 , 'p2',2 , 'p4',4 , 'p10',10 );
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=utRayleigh(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmSigma)
            % Computations specific to the Rayleigh distribution.
            testCase.Dist = Rayleigh(parmSigma);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            
            SetupXs(testCase,41,500);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
%           testCase.MGFMom2RelTol = 0.05;
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % utRayleigh


