classdef utExGamma < utContinuous
    
    properties (ClassSetupParameter)
        % Avoid similar rates for G and E, as this produces nearly singular info matrices in estimation.
        % Seems much slower when parmGscale > parmEmean & can bomb if discrepancy is large.
        parmGshape = struct( 'p1_5',1.5 , 'p120',120, 'p30', 30, 'p10',10);
        parmGscale = struct( 'p400',400 , 'p4',4    , 'p25', 25, 'p75',75);
        parmEmean  = struct( 'p350',350 , 'p200',200, 'p153',153,'p68',68);
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utExGamma(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmGshape,parmGscale,parmEmean)
            % Computations specific to the ExGamma distribution.
            testCase.Dist = ExGamma(parmGshape,parmGscale,parmEmean);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,2000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utExGamma


