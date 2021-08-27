classdef utLaplace < utContinuous
    
    properties (ClassSetupParameter)
        parmLocation = struct( 'n100',-100 , 'p0',0 , 'p5',5   , 'p25',25   , 'p250',250 );
        parmScale    = struct( 'p_5',.5    , 'p1',1 , 'p10',10 , 'p100',100 , 'p11',11   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utLaplace(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmLocation,parmScale)
            % Computations specific to the Laplace distribution.
            testCase.Dist = Laplace(parmLocation,parmScale);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,41,500);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmLocation;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            
            testCase.Expected.Median = parmLocation;
            testCase.Expected.PctileSkew75 = 0;
            testCase.Expected.PctileSkew90 = 0;
            testCase.Expected.PctileSkew99 = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.ParmEstAbsTol = 0.002;
            testCase.ParmEstRelTol = 0.002;
            testCase.RawMomentAbsTol(4) = 0.1;
            testCase.RawMomentRelTol(4) = 0.1;
            testCase.MLParmTolSE = 0.2;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utLaplace


