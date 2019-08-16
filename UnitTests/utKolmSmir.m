classdef utKolmSmir < utContinuous;
    
    properties (ClassSetupParameter)
        parmDF = struct( 'p10',10 , 'p40',40 ,   'p160', 160   ,  'p640', 640, 'p2560', 2560  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utKolmSmir(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmDF)
            % Computations specific to the KolmSmir distribution.
            testCase.Dist = KolmSmir(parmDF);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.0002);
            testCase.MGFMom2AbsTol = 0.001;  % RELAX
            testCase.MGFMom2RelTol = 0.001;
            testCase.MLParmTolSE = nan;  % SE not computed with integer parameters.
 
            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utKolmSmir


