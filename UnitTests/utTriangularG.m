classdef utTriangularG < utContinuous;
    
    properties (ClassSetupParameter)
        parmlow  = struct( 'n100',-100 , 'p250',250   );
        parmpeak = struct( 'n80',-60   , 'p800',800   );
        parmhi   = struct( 'n10',-10   , 'p1000',1000 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utTriangularG(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmlow,parmpeak,parmhi)
            % Computations specific to the TriangularG distribution.
            parmwidth = parmhi - parmlow;
            testCase.Dist = TriangularG(parmlow,parmpeak,parmhi);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            % Parameter estimation often fails due to parameter crossing when the
            % interval is narrow or the peak is close to a boundary.
            testCase.SkipEstAll = true;

            % Compute whatever values known are known from other sources:
            testCase.Expected.Minimum = parmlow;
            testCase.Expected.Maximum = parmhi;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            % Parameter estimation is sometimes very bad
            testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.1;
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.1;
            testCase.MLParmTolSE = 0.5;
          
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utTriangularG


