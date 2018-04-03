classdef utExtrVal2 < utContinuous;
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'p_5',.5 , 'p1',1   , 'p30',30   , 'p10',10 , 'p500',500 );
        parmscale = struct( 'p5',5   , 'p3',3   , 'p2',2     , 'p8',8   , 'p20',20 );
        parmshape = struct( 'p5',5   , 'p10',10 , 'p3',3     , 'p4',4   , 'p11',11 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExtrVal2(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmscale,parmshape)
            % Computations specific to the ExtrVal2 distribution.
            testCase.Dist = ExtrVal2(parmmu,parmscale,parmshape);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
%            testCase.MLParmTolSE = 0.5;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExtrVal2


