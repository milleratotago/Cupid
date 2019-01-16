classdef DemoFitConstrainedFns
    % This class is just a container for various constraint functions
    % used in illustrating constrained fitting.
    
    methods (Static)
        
        function [Dists, Penalty] = Fn0(Dists,Parms)
            % Example: Function to implement the constraint sigma>0.1*mu.
            %
            % Parms is a vector of 2 suggested parameter values from fminsearch.
            % We will treat the first parameter as the suggested mu
            % and the second parameter as the sqrt(suggested sigma).
            % (That is, we will square the second parameter to get the
            % suggested sigma, thereby insuring that sigma is positive.)
            
            % Compute suggested parameter values based on fminsearch's suggestions.
            Mu = Parms(1);
            SD = Parms(2)^2;
            
            % Reset the distribution to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            Dists{1}.ResetParms([Mu SD]);
            
            % Compute a penalty to be added to the error score
            % when the parameters do not satisfy the constraint.
            % Since fminsearch is trying to minimize error, it will
            % tend to search for parameters for which Penalty is zero.
            
            if SD<Mu*0.1
                % Constraint is not satisfied. Penalty increases
                % to the extent that SD is smaller than 0.1*mu.
                Penalty = 100000*(SD-Mu*0.1)^2;
            else
                % Constraint is satisfied, so no penalty.
                Penalty = 0;
            end
            
        end
        
        function [Dists, Penalty] = Fn1(Dists,Parms)
            % Example: Function to implement the constraint that two
            % distributions have the same sigma.
            %
            % Parms: vector of 3 values: suggested values for mu1, mu2, and sqrt(sigma)
            Sigma = Parms(3)^2;
            
            % Reset the distributions to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            Dists{1}.ResetParms([Parms(1) Sigma]);  % Suggested parameters for distribution 1
            Dists{2}.ResetParms([Parms(2) Sigma]);  % Suggested parameters for distribution 2
            
            Penalty = 0;  % No penalty is needed for this case, because it is impossible
            % for fminsearch to suggest values violating the constraint.
        end
        
        function [Dists, Penalty] = Fn2(Dists, Parms)
            % Example constraint function to implement constraints
            % among 3 to-be-fitted distributions. The distributions are
            % Normal, Exponential, and ExGaussian, and the constraint is that
            % the parameters of the ExGaussian equal the corresponding
            % parameters of the normal and exponential.
            %
            % Parms: vector of these 3 values:
            NormalMu = Parms(1);
            NormalSD = Parms(2)^2;  % make sure it is positive
            ExponentialMean = Parms(3)^2;  % make sure it is positive
            
            % Reset the distributions to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            Dists{1}.ResetParms([NormalMu NormalSD]);  % Suggested parameters for distribution 1
            Dists{2}.ResetParms(ExponentialMean);      % Suggested parameters for distribution 2
            Dists{3}.ResetParms([NormalMu NormalSD ExponentialMean]); % ... for distribution 3
            
            Penalty = 0;
            
        end
        
        function [Dists, Penalty] = Fn4(Dists,Parms)
            % Example constraint function for the situation described by Bausenhart et al.
            % The problem is to fit two logistic(location,spread) distributions
            % where there are only 3 free parameters: the constraint allows
            % us to compute the location of the second logistic from the
            % location of the other, plus the two spreads.
            
            % Parms: vector of these 3 values suggested by fminsearch:
            Location1 = Parms(1);
            Spread1 = Parms(2)^2;  % make sure it is positive
            Spread2 = Parms(3)^2;  % make sure it is positive
            
            s = 50;  % an arbitrary constant that is involved in the definition of the constraint.
            Location2 = Spread2/Spread1 * (s - Location1) + s;
            
            % Reset the distributions to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            Dists{1}.ResetParms([Location1 Spread1]);   % Suggested parameters for distribution 1
            Dists{2}.ResetParms([Location2 Spread2]);   % Suggested parameters for distribution 2
            
            Penalty = 0;
            
        end
        
        function [Dists, Penalty] = Fn5(Dists,Parms)
            % Example constraint function to implement the constraint
            % that the estimated distributions should have the same
            % means and SDs.
            %
            % Parms: vector of the parameter values:
            
            % Reset the distributions to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            Dists{1}.ResetParms(Parms(1:2));
            Dists{2}.ResetParms(Parms(3:4));
            
            % Compute a penalty reflecting the differences in means
            % and SDs of the two distributions
            MeanErrSq = (Dists{1}.Mean - Dists{2}.Mean)^2;
            SDErrSq = (Dists{1}.SD - Dists{2}.SD)^2;
            
            Penalty = 100000 * (MeanErrSq + SDErrSq);
            
        end
        
        function [Dists, Penalty] = Fn6(Dists, Parms)
            % Example constraint function to implement constraints
            % among 3 to-be-fitted distributions: The third distribution must
            % be a mixture of the first two.
            %
            % Parms: vector of these values:
            MixP = Parms(1)^2/(1+Parms(1)^2);  % Make sure it is between 0 and 1
            Normal1Mu = Parms(2);
            Normal1SD = Parms(3)^2;  % make sure it is positive
            Normal2Mu = Parms(4);
            Normal2SD = Parms(5)^2;  % make sure it is positive
            
            % Reset the distributions to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            Dists{1}.ResetParms([Normal1Mu Normal1SD]);
            Dists{2}.ResetParms([Normal2Mu Normal2SD]);
            Dists{3}.ResetParms([MixP Normal1Mu Normal1SD Normal2Mu Normal2SD]);
            
            Penalty = 0;
            
        end
        
        function [Dists, Penalty] = Fn7(Dists,Parms)
            % Example constraint function to implement the constraint
            % that the estimated distributions should have the same
            % means and SDs.
            %
            % Parms: vector of the parameter values:
            
            % Reset the distributions to have the parameter values
            % corresponding to fminsearch's suggestions in Parms.
            
            % The input Parms values are any real numbers that fminsearch decides to try.
            % The next two statements convert these arbitrary reals into legal parameter
            % values for the two distributions.
            
            Parms1 = Dists{1}.RealsToParms(Parms(1:Dists{1}.NDistParms));
            Parms2 = Dists{2}.RealsToParms(Parms(Dists{1}.NDistParms+1:end));
            
            Dists{1}.ResetParms(Parms1);
            Dists{2}.ResetParms(Parms2);
            
            % Compute a penalty reflecting the differences in means
            % and SDs of the two distributions
            MeanErrSq = (Dists{1}.Mean - Dists{2}.Mean)^2;
            SDErrSq = (Dists{1}.SD - Dists{2}.SD)^2;
            
            Penalty = 100000 * (MeanErrSq + SDErrSq);
            
        end
        
        
    end % methods (Static)
    
end % classdef
