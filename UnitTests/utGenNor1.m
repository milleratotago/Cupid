classdef utGenNor1 < utContinuous
    
    properties (ClassSetupParameter)
        parmmu    = struct( 'n100',-100 , 'n1',-1    , 'p0',0     , 'p5',5     , 'p25',25   , 'p250',250 );
        parmalpha = struct( 'p_75',.75  , 'p_9',.9   , 'p1',1     , 'p10',10   , 'p100',100 , 'p11',11   );
        %  Beta 1 & 2 should match Laplace and Normal
        parmbeta  = struct( 'p1_0',1.0  , 'p2_0',2.0 , 'p1_3',1.3 , 'p2_4',2.4 , 'p1_5',1.5 , 'p3' ,3   );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utGenNor1(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmmu,parmalpha,parmbeta)
            % Computations specific to the distribution.
            testCase.Dist = GenNor1(parmmu,parmalpha,parmbeta);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            testCase.SkipEstMoment = true;  % Moments do not provide enough info to constrain parameters
            SetupXs(testCase,100,2000);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = parmmu;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            testCase.Expected.Median = parmmu;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.001);
            testCase.RawSkewAbsTol = max(0.001, 0.01 * (abs(parmmu) + parmalpha));  % RawSkewness integral is problematic for large mu and/or alpha
            testCase.RawMomentAbsTol(4) = max(0.001, 0.01 * (abs(parmmu) + parmalpha));  % RawSkewness integral is problematic for large mu and/or alpha
            testCase.ParmEstAbsTol = 0.02;
            testCase.ParmEstRelTol = 0.02;

            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function ComparePDFs(testCase)
            % Check matches to known PDFs
            if testCase.Dist.Beta == 1
                tempU = Laplace(testCase.Dist.Mean,testCase.Dist.Alpha);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            elseif testCase.Dist.Beta == 2
                tempU = Normal(testCase.Dist.Mean,testCase.Dist.SD);
                tempPDF = tempU.PDF(testCase.xvalues);
                testCase.verifyEqual(testCase.Computed.PDF,tempPDF,'AbsTol',1e-6,'Does not match expected PDF values.');
            end
        end
        
    end
    
end  % utGenNor1


