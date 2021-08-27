classdef uttNoncentral < utContinuous
    
    properties (ClassSetupParameter)
        parmdf     = struct( 'p6',6       , 'p12',12 , 'p24',24, 'p120',120);
        parmnoncen = struct( 'p_005',.005 , 'p_3',.3 , 'p1',1  , 'p_1',.1);
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end
    
    methods
        
        function testCase=uttNoncentral(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmdf,parmnoncen)
            % Computations specific to the tNoncentral distribution.
            testCase.Dist = tNoncentral(parmdf,parmnoncen);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            
            testCase.EstParmCodes = 'fr';
            
            SetupXs(testCase,40,1000);
            
            % testCase.Expected.PDF = nctpdf( testCase.xvalues,parmdf,parmnoncen );
            % testCase.Expected.CDF = nctcdf( testCase.xvalues,parmdf,parmnoncen );
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.01);  % Not numerically very accurate
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % uttNoncentral


