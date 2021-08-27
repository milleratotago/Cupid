classdef utProdUni0p < utContinuous
    % Special test class for MultTrans implementation of distribution of product of K Uniform(0,p)
    properties (ClassSetupParameter)
        parmK   = struct( 'K1',1     , 'K2',2     , 'K3',3     , 'K5',5  );  % Big numerical problems already at K = 10
        parmp   = struct( 'p_37',0.37, 'p_38',0.38, 'p_49',0.49, 'p1_5',1.5  );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utProdUni0p(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmK,parmp)
            testCase.Dist = MultTrans(ProdUni01(parmK),parmp^parmK);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,40,200);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);

            utGenericMethodSetup(testCase);   % Initialize many standard computations
        
        end
        
    end  % TestClassSetup

        
end  % utProdUniformK


