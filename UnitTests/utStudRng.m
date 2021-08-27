classdef utStudRng < utContinuous
    
    properties (ClassSetupParameter)
        parmDF = struct( 'p50',50 , 'p70',70 , 'p90',90  );
        parmR  = struct( 'p3',3   , 'p4',4   , 'p5',5  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utStudRng(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmDF,parmR)
            % Computations specific to the t distribution.
            testCase.Dist = StudRng(parmDF,parmR);
            testCase.Dist.UseSplineCDFOn(501);   % In the interest of speed
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            testCase.SkipEstAll = true;  % Skip in interests of speed

            SetupXs(testCase,40,200);
            
            testCase.MGFh = 1.0E-4; % Too large or too small produces numerical errors.

            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Not very accurate
            testCase.RawMomentAbsTol(3) = .07;  % Especially inaccurate.
            testCase.RawMomentAbsTol(4) = .10;
            testCase.CenMomentAbsTol(2) = .03;
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utt


