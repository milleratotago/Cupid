classdef utExtrValGen < utContinuous;
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'p_5',.5 , 'p1',1   , 'p30',30  , 'p10',10  , 'p500',500 );
        parmscale = struct( 'p5',5   , 'p3',3   , 'p_7',.7  , 'p8',8    , 'p1',1     );
        parmshape = struct( 'p_5',.5 , 'p_2',.2 , 'n_1',-.1 , 'n_4',-.4 , 'p0',0 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExtrValGen(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmscale,parmshape)
            % Computations specific to the ExtrValGen distribution.
            testCase.Dist = ExtrValGen(parmmu,parmscale,parmshape);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.EstParmCodes = 'frr';

            SetupXs(testCase,40,500);

            SetTolerances(testCase,0.005);
            testCase.CenMomentAbsTol(2) = 0.01;
            testCase.RawMomentAbsTol(3) = 0.01;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExtrValGen


