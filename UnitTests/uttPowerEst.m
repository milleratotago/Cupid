classdef uttPowerEst < utContinuous;
    
    properties (ClassSetupParameter)
        % Parm values to be combined sequentially.
        parmCase = {1 2 3 4};
    end
    
    properties
        Dummy1, Dummy2  % Provide the dummy variables that were used to make parent classes abstract.
        ThisCase
    end
    
    methods
        
        function testCase=uttPowerEst(varargin)  % Constructor
            testCase=testCase@utContinuous(varargin{:});
        end
        
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        
        function ClassSetup(testCase,parmCase)
            % Computations specific to the tPowerEst distribution.
            testCase.ThisCase  = parmCase;
            switch parmCase
                                  % tPowerEst(TrueDelta,sigma,alpha,n1,[n2])
                case 1
                    testCase.Dist = tPowerEst(.4,1,.05,20);
                case 2
                    testCase.Dist = tPowerEst(.2,1,.05,50);
                case 3
                    testCase.Dist = tPowerEst(.4,1,.05,20,20);
                case 4
                    testCase.Dist = tPowerEst(.2,1,.05,50,50);
            end
            fprintf('\nInitialized %s\n',testCase.Dist.StringName)
            %             testCase.Dist.SearchOptions.MaxFunEvals = 30000;
            %             testCase.Dist.SearchOptions.MaxIter = 20000;
            
            testCase.SkipAllEst = true;  % Parameter estimation is too slow
            testCase.Dist.UseSplineCDFOn(200);
            testCase.Dist.UseSplinePDFOn(200);

            SetupXs(testCase,40,1000);
            
            % Adjust tolerances as appropriate for this distribution & parameters:
            SetTolerances(testCase,0.005);
            % if testCase.ThisCase==1
            %     testCase.ParmEstRelTol(:)=.05;
            % end
            % if testCase.ThisCase==2
            %     testCase.ParmEstRelTol(:)=.2;
            %     testCase.ParmEstAbsTol(:)=.2;
            % end
            % if testCase.ThisCase==3
            %     testCase.ParmEstRelTol(:)=.1;
            %     testCase.ParmEstAbsTol(:)=.1;
            % end
            
            utGenericMethodSetup(testCase);   % Initialize many standard computations
            
        end
        
    end  % TestClassSetup
    
end  % uttPowerEst


