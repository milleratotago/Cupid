classdef FNoncentral < dContinuous
    % FNoncentral(DfNum,DfDenom,Noncentrality)
    
    properties(SetAccess = protected)
        dfnum, dfdenom, noncen,
        Sqrtnoncen, dfSearchMax, noncenMax, Tolerance, CentralF,
        UseCentral;
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.GT2Real(1,Parms(2)) NumTrans.GT2Real(0,Parms(3))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2GT(1,Reals(2)) NumTrans.Real2GT(0,Reals(3))];
        end
        
    end
    
    methods
        
        function obj=FNoncentral(varargin)
            obj=obj@dContinuous('FNoncentral');
            obj.ParmTypes = 'iir';
            obj.DefaultParmCodes = 'iir';
            obj.NDistParms = 3;
            obj.CentralF = F;
            obj.CDFNearlyZero = 1e-8;
            obj.CDFNearlyOne = 1 - obj.CDFNearlyZero;
            %            obj.ComputeBounds = true;
            %            obj.IntMinSteps = 6;
            obj.Tolerance = 1.0E-8;  % determines CDF convergence.
            obj.noncenMax = 2500;    % A bound to keep noncen from running away in searches.
            obj.dfSearchMax = 5000;  % A bound to keep dfs from running away in searches.
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('FNoncentral:Constructor', ...
                        'FNoncentral constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.dfnum = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.dfdenom = VerifyIntegerGE(obj,1,newparmvalues(2));
            obj.noncen = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newdfnum   = ifelse(ParmCodes(1)=='f', obj.dfnum, obj.dfnum+1);
            newdfdenom = ifelse(ParmCodes(2)=='f', obj.dfdenom, obj.dfdenom+1);
            newnoncen  = ifelse(ParmCodes(3)=='f', obj.noncen,1.1*obj.noncen);
            obj.ResetParms([newdfnum newdfdenom newnoncen]);
        end
        
        function []=ReInit(obj)
            assert(obj.noncen>0,'FNoncentral noncen must be > 0.');
            obj.Initialized = true;
            obj.UseCentral = obj.noncen < 1.0E-10; % Brown's criterion for using central F.
            if obj.UseCentral
                obj.CentralF.ResetParms([obj.dfnum obj.dfdenom]);
                obj.LowerBound = obj.CentralF.LowerBound;
                obj.UpperBound = obj.CentralF.UpperBound;
            else
                % Need better setting of bounds.
                obj.LowerBound = obj.XNearlyZero;
                obj.UpperBound = 0.5;
                Thiscdf = -0.5;
                Oldcdf = -1;
                while (Thiscdf < obj.CDFNearlyOne) && (Thiscdf > Oldcdf)
                    obj.UpperBound = obj.UpperBound * 2;
                    Oldcdf = Thiscdf;
                    Thiscdf = obj.CDF(obj.UpperBound-obj.XNearlyZero);
                end
                % fprintf('A: %6.3f %10.3f\n',obj.LowerBound,obj.UpperBound);
                if Thiscdf < obj.CDFNearlyOne
                    warning(['CDF of FNoncentral asymptote is below CDFNearlyOne at ' num2str(Thiscdf)]);
                    obj.UpperBound = obj.UpperBound / 2;
                    while obj.CDF(obj.UpperBound) >= Thiscdf
                        obj.UpperBound = 0.99 * obj.UpperBound;
                    end
                    obj.UpperBound = obj.UpperBound / 0.99;
                else
                    obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
                end
                % fprintf('B: %6.3f %10.3f\n',obj.LowerBound,obj.UpperBound);
                if obj.CDF(obj.LowerBound) < obj.CDFNearlyZero
                    obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
                end
                % fprintf('C: %6.3f %10.3f\n',obj.LowerBound,obj.UpperBound);
            end
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)   % Wikipedia has it, but so does MATLAB
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = ncfpdf(X(InBounds),obj.dfnum,obj.dfdenom,obj.noncen);
        end

        function thiscdf=CDF(obj,X)
            % Based on a Pascal translation of Barry W. Brown's cumfnc.f
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            if obj.UseCentral
                % Handle case in which the non-centrality parameter is (essentially) 0.
                thiscdf = obj.CentralF.CDF(X);
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    % thiscdf(i) = 1 - betainc(obj.dfdenom/(obj.dfdenom+obj.dfnum*X(i)),obj.dfdenom/2,obj.dfnum/2);
                    xnonc = obj.noncen/2;
                    % Calculate the central term of the poisson weighting factor
                    icent = floor(xnonc+0.5);   % Rounding assumed here.
                    if icent == 0
                        icent = 1;
                    end
                    % Compute central weight term
                    centwt = exp(-xnonc+icent*log(xnonc)-gammaln(icent+1));
                    % Compute central incomplete beta term.  Assure that minimum of arg to beta and 1 - arg is computed accurately.
                    ddfn = obj.dfnum;
                    ddfd = obj.dfdenom;
                    df = X(iel);
                    prod = ddfn * df;
                    dsum = ddfd + prod;
                    yy = ddfd / dsum;
                    if yy > 0.5
                        xx = prod / dsum;
                        yy = 1 - xx;
                    else
                        xx = 1 - yy;
                    end
                    
                    dbetdn = betainc(xx,ddfn*0.5+icent,ddfd*0.5);
                    
                    betdn = dbetdn;
                    arg = xx;
                    omarg = yy;
                    adn = obj.dfnum/2 + icent;
                    aup = adn;
                    b = obj.dfdenom/2;
                    betup = betdn;
                    sum = centwt * betdn;
                    
                    % Now sum terms backward from icent until convergence or all done
                    
                    xmult = centwt;
                    i = icent;
                    jtemp = gammaln(adn+b);
                    jtemp = jtemp-gammaln(adn+1.0);
                    jtemp = jtemp-gammaln(b);
                    dnterm = exp(jtemp+adn*log(arg)+b*log(omarg));
                    
                    while (~QSmall(xmult*betdn)) || (i > 0)
                        xmult = xmult * (i/xnonc);
                        i = i - 1;
                        adn = adn - 1;
                        dnterm = (adn+1)/ ((adn+b)*arg)*dnterm;
                        betdn = betdn + dnterm;
                        sum = sum + xmult*betdn;
                    end
                    
                    i = icent + 1;
                    
                    % Now sum forward until convergence
                    
                    xmult = centwt;
                    if aup-1+b == 0
                        jtemp = gammaln(aup);
                        jtemp = jtemp + gammaln(b);
                        upterm = exp( -jtemp + (aup-1)*log(arg) + b*log(omarg) );
                    else
                        jtemp = gammaln(aup-1+b);
                        jtemp = jtemp - gammaln(aup);
                        jtemp = jtemp - gammaln(b);
                        upterm = exp( jtemp + (aup-1)*log(arg) + b*log(omarg) );
                    end
                    
                    Converged = false;
                    while ~Converged
                        xmult = xmult * (xnonc/i);
                        i = i + 1;
                        aup = aup + 1;
                        upterm = (aup+b-2) * arg / (aup-1) * upterm;
                        betup = betup - upterm;
                        sum = sum + xmult * betup;
                        Converged = QSmall(xmult*betup);
                    end
                    thiscdf(iel) = sum;
                end
            end
            function thiscdf=QSmall(X)
                thiscdf = (sum < 1.0E-20) || (X < obj.Tolerance*sum);
            end
        end
        
        function thisval=Mean(obj) % JKB 2nd ed vol 2 p 481
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.dfdenom > 2
                thisval = obj.dfdenom * (obj.dfnum + obj.noncen) / (obj.dfnum * (obj.dfdenom-2) );
            else
                thisval = inf;
            end
        end
        
        function thisval=Variance(obj)   % JKB 2nd ed vol 2 p 481:
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if obj.dfdenom <= 4
                thisval = inf;  %  'Requested infinite variance of FObj.Noncentral distribution.'
            else
                thisval = 2 * (obj.dfdenom/obj.dfnum)^2 * ( (obj.dfnum+obj.noncen)^2 + (obj.dfnum+2*obj.noncen)*(obj.dfdenom-2) ) / ( (obj.dfdenom-2)^2*(obj.dfdenom-4) );
            end
        end
        
    end  % methods
    
end  % class FNoncentral




