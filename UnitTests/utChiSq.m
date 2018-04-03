classdef utChiSq < utContinuous;
    
    properties (ClassSetupParameter)
        parmDF = struct( 'p1',1 , 'p2',2 , 'p4',4 ,  'p8', 8   ,  'p16', 16   ,  'p32', 32  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utChiSq(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmDF)
            % Computations specific to the ChiSq distribution.
            testCase.Dist = ChiSq(parmDF);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);

            % Adjust tolerances as appropriate for this distribution & parameters:
            if parmDF > 1
                SetTolerances(testCase,0.0002);
            else
                SetTolerances(testCase,0.002);
                testCase.CenMomentRelTol(3) = .005;
                testCase.RawMomentRelTol(4) = 0.01;
                testCase.KurtRelTol = 0.05;
                testCase.RawIntAbsTol(3) = 0.01;
                testCase.RawIntRelTol(4) = 0.015;
                testCase.RawIntRelTol(5) = 0.04;
            end
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters.
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utChiSq


