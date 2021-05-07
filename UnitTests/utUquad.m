classdef utUquad < utContinuous;
    
    properties (ClassSetupParameter)
        parmlow = struct( 'n100',-100 , 'n1',-1   , 'p0',0 , 'p5',5   , 'p250',250   );
        parmhi  = struct( 'n10',-10   , 'p_9',.9  , 'p1',1 , 'p10',10 , 'p1000',1000 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utUquad(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmlow,parmhi)
            % Computations specific to the Uquad distribution.
            testCase.Dist = Uquad(parmlow,parmhi);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,101,10001);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = (parmlow+parmhi)/2;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            
            testCase.Expected.Median = testCase.Expected.Mean;
            testCase.Expected.Minimum = parmlow;
            testCase.Expected.Maximum = parmhi;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
            testCase.RawMomentAbsTol(4) = .004;
            testCase.MLParmTolSE = 1.5;   % ML parameter estimation is not great
            % Unfortunately, moment accuracy is poorer (so tolerances must be larger)
            %  when the width is very small or large and when the mean is larger.
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utUquad


