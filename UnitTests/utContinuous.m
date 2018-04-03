classdef utContinuous < utGeneric;
    
    % An abstract class for checking properties of continuous distributions:
    %    checking pdf/cdf values against known values, and
    %    internal consistency tests that should be passed by _every_ continuous probability distribution.
    
    properties (Abstract)
        % Use something Abstract so that matlab.unittest.TestSuite.fromFolder
        % will NOT attempt to define tests from this file.
        Dummy2
    end
    
    methods
        
        function testCase=utContinuous(varargin)  % Constructor
            testCase=testCase@utGeneric(varargin{:});
            % utGeneric initializes these properties:
%             testCase.HighestMoment = 4;
%             testCase.ChiSqNRands = ;
%             testCase.ChiSqNBins = ;
%             testCase.ChiSqNTries = ;
%             testCase.ChiSqCriticalp = ;
            % utContinuous initializes these properties:
        end

        function SetupXs(testCase,nxs,nxmles)

            % Set up some X values at which PDF, CDF, etc should be evaluated (also used in non-MLE estimation);
            testCase.xvalues = testCase.Dist.InverseCDF( (1:2:(2*nxs-1)) / (2*nxs) );

            % Set up some X values for which MLE should return (very close to) the true parameters:
            testCase.xMLE = testCase.Dist.InverseCDF( (1:2:(2*nxmles-1)) / (2*nxmles) );

        end
            
    end  % regular methods
    
    methods (Test)
        
        function PDFIntegralPieces(testCase)
            
            % Check for consistency of CDF with integral of PDF:
            NSegments = numel(testCase.xvalues);
            IntegralPieces = zeros(1,NSegments);
            LowerLimit = testCase.Dist.LowerBound;
            for iSegment=1:NSegments
                UpperLimit = testCase.xvalues(iSegment);
                IntegralPieces(iSegment) = testCase.Dist.IntegralXToNxPDF(LowerLimit,UpperLimit,0);
                LowerLimit = UpperLimit;
            end
            IntegralPieces = cumsum(IntegralPieces);
            testCase.verifyEqual(testCase.Computed.CDF,IntegralPieces,'AbsTol',testCase.CDFAbsTol,'RelTol',testCase.CDFRelTol,'CDF values are not consistent with integrated PDF');
            
        end  % PDFIntegralPieces
        
        
    end  %     methods (Test)
    
end

