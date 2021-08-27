classdef utDblMon < utContinuous
    
    properties (ClassSetupParameter)
        parm1tzero   = struct( 'p150',150 , 'p200',200 ,   'p500',500   );
        parm2delta   = struct(  'p45',45  ,  'p21',21  ,    'p80',80    );
        parm3epsilon = struct(  'p15',15  ,  'p12',12  , 'p100_1',100.1 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utDblMon(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parm1tzero,parm2delta,parm3epsilon)
            % Computations specific to the DblMon distribution.
            testCase.Dist = DblMon(parm1tzero,parm2delta,parm3epsilon);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,41,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);
            testCase.RawMomentRelTol(2) = 0.05;
            testCase.RawMomentRelTol(3) = 0.10;
            testCase.RawMomentRelTol(4) = 0.15;
            testCase.ParmEstRelTol = 0.1;
            testCase.MLParmTolSE = 1.5;
%             testCase.RawMomentRelTol(3) = 0.01;
            testCase.KurtRelTol = 0.02;
           
%             testCase.MLParmTolSE = 0.05;   % ML parameter estimation is bad
%             testCase.MGFMom2RelTol = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utDblMon


