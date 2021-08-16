classdef utConditXLTY < utContinuous
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4 5};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=utConditXLTY(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the ConditXLTY distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                case 1
                    testCase.Dist = ConditXLTY(Triangular(1,2),Uniform(0.5,1.5));
                    testCase.EstParmCodes = 'fffr';
%                    testCase.SkipEstAll = true;
                case 2
                    testCase.Dist = ConditXLTY(Beta(18,5),Beta(21,2));
                    testCase.EstParmCodes = 'fffr';
%                    testCase.SkipEstAll = true;
                case 3
                    testCase.Dist = ConditXLTY(RNGamma(10,.1),Uniform(1,200));
                    testCase.EstParmCodes = 'frff';
                case 4  % Check known mean from Rolf (see fn CompareMeans)
                    testCase.Dist = ConditXLTY(Exponential(0.01),AddTrans(Exponential(0.025),-40));
                    testCase.EstParmCodes = 'frf';
                case 5  % Check known mean from Rolf (see fn CompareMeans)
                    testCase.Dist = ConditXLTY(Exponential(0.01),AddTrans(Exponential(0.025),40));
                    testCase.EstParmCodes = 'frf';
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            %             testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            %             testCase.Dist.SearchOptions.MaxIter = 20000;
            
            % testCase.Dist.UseSplineCDFOn(200);
            % testCase.Dist.UseSplinePDFOn(200);

            SetupXs(testCase,40,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            if testCase.ThisCase==1
%                testCase.ParmEstRelTol(:)=.05;
            end
            if testCase.ThisCase==2
%                testCase.ParmEstRelTol(:)=.2;
%                testCase.ParmEstAbsTol(:)=.2;
            end
            if testCase.ThisCase==3
%                testCase.ParmEstRelTol(:)=.1;
%                testCase.ParmEstAbsTol(:)=.1;
            end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
    methods (Test)
        
        function CompareMeans(testCase)
            % Check matches to known means.
            % These known means are from Rolf Ulrich's habilitation (1988) p 116, eqn 2.34
            % "Mathematisierte Theorienbildung in der kognitiven Psychologie". Tübingen: Habilitationsschrift
            if testCase.ThisCase == 4
                % X is exponential, Y is shifted exponential with any negative shift d
                % Note conditional mean of X does not depend on d!!!
                lambda_x = testCase.Dist.BasisRV1.rate;
                lambda_y = testCase.Dist.BasisRV2.BasisRV.rate;
                tempMean = 1 / (lambda_x + lambda_y);
                testCase.verifyEqual(testCase.Computed.Mean,tempMean,'AbsTol',1e-4,'Does not match expected mean.');
            elseif testCase.ThisCase == 5
                % X is exponential, Y is shifted exponential with positive shift d
                lambda_x = testCase.Dist.BasisRV1.rate;
                lambda_y = testCase.Dist.BasisRV2.BasisRV.rate;
                d = testCase.Dist.BasisRV2.Addend;
                tempMeanNumerator = 1 / lambda_x - ...
                      lambda_y*exp(-lambda_x*d) / (lambda_x + lambda_y)^2 * ( lambda_y / lambda_x + 2 + d*(lambda_x + lambda_y) );
                tempMeanDenominator = 1 - exp(-lambda_x*d) * lambda_y / (lambda_x + lambda_y);
                tempMean = tempMeanNumerator / tempMeanDenominator;
                testCase.verifyEqual(testCase.Computed.Mean,tempMean,'AbsTol',1e-4,'Does not match expected mean.');
            end
        end
        
    end
    
end  % utConditXLTY

