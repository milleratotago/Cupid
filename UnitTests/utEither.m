classdef utEither <  utContinuous & utDiscrete  % NWJEFF: Never used. Could combine eg utAddTrans and utAddTransDisc
    
    methods
        
        function testCase=utEither(varargin)  % Constructor
            testCase@utContinuous(varargin{:});
            testCase@utDiscrete(varargin{:});
        end

    end  % regular methods
    
    methods (Test)
        
        function PDFIntegralPieces(testCase)
            switch testCase.Dist.DistType
                case 'c'
                    PDFIntegralPieces@utContinuous(testCase);
                case 'd'
                    PDFIntegralPieces@utDiscrete(testCase);
            end
        end
        
    end
    
end  % utEither


