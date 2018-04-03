classdef utCosine < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmmu        = struct( 'n10',-10  , 'p2',2 , 'p27',27 , 'p100',100  );
        parmhalfwidth = struct( 'p_05',.05 , 'p2',2 , 'p1',1   , 'p12',12    );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utCosine(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmhalfwidth)
            % Computations specific to the Cosine distribution.
            testCase.Dist = Cosine(parmmu,parmhalfwidth);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            testCase.MGFh = 1.0E-4*parmhalfwidth; % Too large or too small produces numerical errors. This value works quite well for normal.

            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmmu;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            
            testCase.Expected.Median = parmmu;
            testCase.Expected.Minimum = parmmu-parmhalfwidth;
            testCase.Expected.Maximum = parmmu+parmhalfwidth;
            testCase.Expected.PctileSkew75 = 0;
            testCase.Expected.PctileSkew90 = 0;
            testCase.Expected.PctileSkew99 = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.MLParmTolSE = 0.3;
            testCase.MGFXRelTol = 0.005;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utCosine


