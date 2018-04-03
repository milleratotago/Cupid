classdef utBinomial < utGeneric;
    
    properties (ClassSetupParameter)
        parmN = struct( 'p2',2   , 'p4',4     , 'p20',20 , 'p100',100 , 'p500',500 );
        parmp = struct( 'p_1',.1 , 'p_83',.83 , 'p_8',.8 , 'p_55',.55 , 'p_73',.73 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utBinomial(varargin)  % Constructor
            testCase=testCase@utGeneric(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmN,parmp)
            % Computations specific to the normal distribution.
            testCase.Dist = Binomial(parmN,parmp);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.xvalues = testCase.Dist.LowerBound:testCase.Dist.UpperBound;
            
            % Set up some X values for which MLE should return the true parameters:
            npoints = 500;
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*npoints-1)) / (2*npoints) );
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmN*parmp;
            testCase.Expected.Variance = parmN*parmp*(1-parmp);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utBinomial


