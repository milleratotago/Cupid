classdef utTriangularGCWP < utContinuous;
    
    properties (ClassSetupParameter)
        parm1center   = struct( 'n100',-100 , 'n1',-1     , 'p0',0     , 'p5',5   , 'p250',250 );
        parm2width    = struct( 'p80',80    , 'p_96',.96  , 'p600',600 , 'p8',8   , 'p90',90   );
        parm3peakprop = struct( 'p_10',.10  , 'p_9',.9    , 'p_5',.5   , 'p_3',.3 , 'p_77',.77 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utTriangularGCWP(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1center,parm2width,parm3peakprop)
            % Computations specific to the TriangularGCWP distribution.
            parmwidth = parm3peakprop - parm1center;
            testCase.Dist = TriangularGCWP(parm1center,parm2width,parm3peakprop);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,500);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            % Parameter estimation is sometimes very bad
            testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.1;
            testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.1;
            testCase.MLParmTolSE = 0.5;
          
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utTriangularGCWP


