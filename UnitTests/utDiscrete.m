classdef utDiscrete < utGeneric;
    
    % An abstract class for checking properties of discrete distributions:
    %    checking pdf/cdf values against known values, and
    %    internal consistency tests that should be passed by _every_ discrete probability distribution.
    
    properties (Abstract)
        % Use something Abstract so that matlab.unittest.TestSuite.fromFolder
        % will NOT attempt to define tests from this file.
        Dummy2
    end
    
    methods
        
        function testCase=utDiscrete(varargin)  % Constructor
            testCase=testCase@utGeneric(varargin{:});
            % utGeneric initializes these properties:
%             testCase.HighestMoment = 4;
%             testCase.ChiSqNRands = ;
%             testCase.ChiSqNBins = ;
%             testCase.ChiSqNTries = ;
%             testCase.ChiSqCriticalp = ;
            % utDiscrete initializes these properties:
        end

        function SetupXs(testCase,~,~)

            % Set up some X values at which PDF, CDF, etc should be evaluated (also used in non-MLE estimation);
            testCase.xvalues = testCase.Dist.XsToPlot;

            % Set up some X values for which MLE should return (very close to) the true parameters:
            testCase.xMLE = testCase.xvalues;

        end
            
    end  % regular methods
    
    methods (Test)
        
        function PDFIntegralPieces(testCase)
            
            % Check for consistency of CDF with sum of PDF:
            IntegralPieces = testCase.Dist.PDF(testCase.xvalues);
            testCase.Computed.CDF = testCase.Dist.CDF(testCase.xvalues);
            IntegralPieces = cumsum(IntegralPieces);
            testCase.verifyEqual(testCase.Computed.CDF,IntegralPieces,'AbsTol',testCase.CDFAbsTol,'RelTol',testCase.CDFRelTol,'CDF values are not consistent with summed PDF');
            
        end  % PDFIntegralPieces
        
        
    end  %     methods (Test)
    
end

