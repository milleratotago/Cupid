classdef utTriangular < utContinuous
    
    properties (ClassSetupParameter)
        parmlow = struct( 'n100',-100 , 'n1',-1   , 'p0',0 , 'p5',5   , 'p250',250   );
        parmhi  = struct( 'n10',-10   , 'n_9',-.9 , 'p1',1 , 'p10',10 , 'p1000',1000 );
    end
    
    properties
       Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
    end

    methods
        
        function testCase=utTriangular(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmlow,parmhi)
            % Computations specific to the Triangular distribution.
            parmwidth = parmhi - parmlow;
            testCase.Dist = Triangular(parmlow,parmhi);
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)

            SetupXs(testCase,101,10001);
            
            % Compute whatever values known are known from other sources:
            testCase.Expected.Mean = (parmlow+parmhi)/2;
            testCase.Expected.RawSkewness = 0;
            testCase.Expected.RelSkewness = 0;
            
            testCase.Expected.Median = testCase.Expected.Mean;
            testCase.Expected.Minimum = parmlow;
            testCase.Expected.Maximum = parmhi;
            testCase.Expected.PctileSkew75 = 0;
            testCase.Expected.PctileSkew90 = 0;
            testCase.Expected.PctileSkew99 = 0;
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.002);
            testCase.RawMomentAbsTol(4) = .004;
            testCase.MLParmTolSE = 1.5;   % ML parameter estimation is not great
            % Unfortunately, moment accuracy is poorer (so tolerances must be larger)
            %  when the width is very small or large and when the mean is larger.
%             RangeFactor = min(parmwidth,1/parmwidth);
%             if RangeFactor>0.5
%                 SkewTol = 0.001;
%             elseif RangeFactor>0.3
%                 SkewTol = 0.002;
%             elseif RangeFactor>0.1
%                 SkewTol = 0.005;
%             else
%                 SkewTol = max(0.03,0.000001*testCase.Expected.Mean/RangeFactor);
%             end
%             testCase.RawSkewAbsTol = SkewTol;
%             testCase.RawSkewRelTol = SkewTol;
%             testCase.ParmEstAbsTol = ones(1,testCase.Dist.NDistParms) * 0.005;
%             testCase.ParmEstRelTol = ones(1,testCase.Dist.NDistParms) * 0.005;
%             testCase.MLErrTol = 0.5;  % Starting error must be greater than EndingError-ThisTolerance.
           
            utGenericMethodSetup(testCase);   % Initialize many standard computations

        end
        
    end  % TestClassSetup
    
end  % utTriangular


