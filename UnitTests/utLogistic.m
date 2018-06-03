classdef utLogistic < utContinuous;
    
    properties (ClassSetupParameter)
        parmmu   = struct( 'n100',-100 , 'n1',-1  , 'p0',0 , 'p5',5 , 'p25',25 , 'p250',250 );
        parmbeta = struct( 'p_05',.05  , 'p_1',.1 , 'p1',1 , 'p5',5 , 'p10',10 , 'p25',25 );  % Roughly half of the SD
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utLogistic(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmbeta)
            % Computations specific to the Logistic distribution.
            testCase.Dist = Logistic(parmmu,parmbeta);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
 
            SetupXs(testCase,40,2000);
            
            testCase.Expected.Median = parmmu;
            testCase.Expected.PctileSkew75 = 0;
            testCase.Expected.PctileSkew90 = 0;
            testCase.Expected.PctileSkew99 = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.0002);
            testCase.RawMomentAbsTol(4) = max( [0.005, 0.001*abs(testCase.Dist.Mean), 0.001*abs(testCase.Dist.SD) ] );
%             testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.001;
%             testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.005;
            testCase.MLParmTolSE = 0.25;   % ML parameter estimation is not great

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utLogistic


