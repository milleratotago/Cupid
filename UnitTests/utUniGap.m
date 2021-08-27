classdef utUniGap < utContinuous
    
    properties (ClassSetupParameter)
        parmt = struct( 'p1',1 , 'p5',5   , 'p10',10 , 'p25',25   , 'p250',250   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utUniGap(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmt)
            % Computations specific to the UniGap distribution.
            testCase.Dist = UniGap(parmt);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,100,2000);
            testCase.xvalues = [testCase.Dist.LowerBound testCase.xvalues testCase.Dist.UpperBound];
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
            testCase.RawMomentAbsTol(4) = 0.01*parmt;  % RawSkewness integral is very problematic
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utUniGap


