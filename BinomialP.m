classdef BinomialP < Binomial
    % BinomialP(N,P) distribution of proportion successes with parameters N_of_trials, Pr_of_success.
    
    methods
        
        function obj=BinomialP(varargin)
            obj=obj@Binomial(varargin{:});
            obj.FamilyName = 'BinomialP';
            BuildMyName(obj);
        end
        
        function []=ReInit(obj)
            obj.ReInit@Binomial;
            obj.DiscreteX = obj.DiscreteX / obj.N;
            obj.DiscreteXmin = obj.DiscreteX - 0.1 / obj.N;
            obj.DiscreteXmax = obj.DiscreteX + 0.1 / obj.N;
            obj.LowerBound = obj.LowerBound / obj.N;
            obj.UpperBound = obj.UpperBound / obj.N;
        end
        
        function thisval=Mean(obj)
            thisval = obj.Mean@Binomial/obj.N;
        end
        
        function thisval=Variance(obj)
            thisval = obj.Variance@Binomial/obj.N^2;
        end
        
        function thisval=MGF(obj,Theta)
            thisval = obj.MGF@Binomial(Theta / obj.N);  % NWJEFF: Also applies to MultTrans & LinearTrans
        end

        function thisval=Random(obj,varargin)
            thisval = obj.Random@Binomial(varargin{:})/obj.N;
        end
        
        function s=EstML(obj,Observations,varargin)
            s = obj.EstML@Binomial(Observations*obj.N,varargin{:});
        end
        
    end  % methods
    
end  % class BinomialP


