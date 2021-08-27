classdef utGenNor2 < utContinuous
    
    properties (ClassSetupParameter)
        parmxi    = struct( 'n100',-100 , 'n1',-1    , 'p0',0     , 'p5',5      , 'p25',25   , 'p250',250 );
        parmalpha = struct( 'p_75',.75  , 'p_9',.9   , 'p1',1     , 'p10',10    , 'p100',100 , 'p11',11   );
        %  Kappa 0 should match Normal
        parmkappa  = struct( 'p0',0     , 'p_2',.2 , 'p_1',0.1 , 'n_4',-0.4 , 'n_5',-0.5 , 'p_3' ,.3   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utGenNor2(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmxi,parmalpha,parmkappa)
            % Computations specific to the distribution.
            testCase.Dist = GenNor2(parmxi,parmalpha,parmkappa);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.SkipEstMoment = true;  % Moments do not provide enough info to constrain parameters
            SetupXs(testCase,100,500);
            
            % Compute whatever values known are known from other sources:
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.RawSkewAbsTol = max(0.001, 0.01 * (abs(parmxi) + parmalpha));  % RawSkewness integral is problematic for large xi and/or alpha
            testCase.RawMomentAbsTol(4) = max(0.001, 0.01 * (abs(parmxi) + parmalpha));  % RawSkewness integral is problematic for large xi and/or alpha
            testCase.ParmEstAbsTol = 0.02;
            testCase.ParmEstRelTol = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.Dist.Kappa == 0
                tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utGenNor2


