classdef NumTrans < handle
    
    % Handy functions for transforming numbers back and forth between scales,
    % often useful when searching for constrained parameters (eg with fminsearch).
    % "real" refers to any real number -inf to +inf.
    % "trans" refers to a number on a constrained scale.
    
    properties (Constant)
        Slope = 0.005
    end % properties (Constant)
    
    methods (Static)
        
        
        % Transform to a number more than some minimum:
        
        function trans = Real2GT(minimum,real)
            % Convert an arbitrary real number to a number greater than some minimum.
            trans = real.^2 + minimum;
        end
        
        function real = GT2Real(minimum,trans)
            % Inverse of Real2GT: Convert a real number with some minimum to an arbitrary real number.
            real = (trans-minimum).^(.5);
        end
        
        
        % Transform to a number less than some maximum:
        
        function trans = Real2LT(maximum,real)
            % Convert an arbitrary real number to a number less than some maximum.
            trans = maximum - real.^2;
        end
        
        function real = LT2Real(maximum,trans)
            % Inverse of Real2LT: Convert a real number with some maximum to an arbitrary real number.
            real = (maximum - trans).^(.5);
        end
        
        
        % Transform to a number within certain bounds:
        
        function trans = Real2Bounded(minimum,maximum,real)
            % Convert an arbitrary real number to a number within certain bounds.
            trans = (maximum - minimum) ./ (exp(-real*NumTrans.Slope) + 1) + minimum;
            if isnan(trans)
                fprintf('Real2Bounded error at min = %f, max = %f, real = %f\n',minimum,maximum,real);
            end
        end
        
        function real = Bounded2Real(minimum,maximum,trans)
            % Inverse of Real2Bounded: Convert a real number within certain bounds to an arbitrary real number.
            real = log( (trans-minimum) ./ (maximum-trans)) / NumTrans.Slope;
            if isnan(real)
                fprintf('Bounded2Real error at min = %f, max = %f, trans = %f\n',minimum,maximum,trans);
            end
        end
        
        
    end  % methods (Static)
    
    
end  % NumTrans
