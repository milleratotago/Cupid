classdef StudRng < dContinuous
    % StudRng(v,r) distribution of the Studentized range statistic used in Tukey's tests.
    % v = degrees of freedom
    % r = number of groups
    
    properties(SetAccess = protected)
        df, r,
        mindf, minr
    end
    
    methods (Static)
        
        function prtrng=cdfTukey(q,v,r)  % NWJEFF: Vectorize?
            %{
Retrieved from https://au.mathworks.com/matlabcentral/fileexchange/37450 on 9 April 2018.
Copyright (c) 2012, Peter Nagy
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
            %}
            % The function returns the cumulative distribution function of the
            % studentized range distribution used in Tukey's HSD test
            % q - observed value of Tukey's test
            % v - degrees of freedom (total number of elements minus number of groups
            % r - number of samples
            pcutj=0.00003;
            pcutk=0.0001;
            step=0.1125; % 0.45 in the original Fortran file
            vmax=120;
            cv1=0.193064705;
            cv2=0.293525326;
            cvmax=0.39894228;
            cv=[0.318309886,-0.268132716e-2, 0.347222222e-2, 0.833333333e-1];
            jmin=12; % 3 in the original Fortran file
            jmax=60; % 15 in the original Fortran file
            kmin=28; % 7 in the original Fortran file
            kmax=60; % 15 in the original Fortran file
            vw=zeros(121); % 30 in the original Fortran file
            qw=zeros(121); % 30 in the original Fortran file
            
            % Minimum and maximum number of steps are controlled by
            % jmin, jmax, kmin and kmax.  Accuracy can be increased
            % by use of a finer grid - Increase sizes of arrays vw
            % and qw, and jmin, jmax, kmin, kmax and 1/step proportionally.
            
            prtrng=0;
            if not(q<=0 || v<1 || r<2)
                g=step*r^-0.2;
                gmid=0.5*log(r);
                r1=r-1;
                c=log(r*g*cvmax);
                if v<=vmax
                    h=step*v^-0.2;
                    v2=v*0.5;
                    switch v
                        case 1
                            c=cv1;
                        case 2
                            c=cv2;
                        otherwise
                            c=sqrt(v2)*cv(1)/(1+((cv(2)/v2+cv(3))/v2+cv(4))/v2);
                    end
                    c=log(c*r*g*h);
                end
                gstep=g;
                qw(1)=-1;
                qw(jmax+1)=-1;
                pk1=1;
                pk2=1;
                k=1;
                while k<=kmax % loop 28
                    gstep=gstep-g;
                    notDone21=true;
                    while notDone21 % loop 21
                        gstep=-gstep;
                        gk=gmid+gstep;
                        pk=0;
                        if not(pk2<=pcutk && k>kmin) % test 26
                            w0=c-gk*gk*0.5;
                            pz=1-normcdf(gk);
                            x=1-normcdf(gk-q)-pz;
                            if x>0, pk=exp(w0+r1*log(x));end
                            if not(v>vmax) % test 26 again
                                jump=-jmax;
                                notDone22=true;
                                while notDone22 % loop 22
                                    jump=jump+jmax;
                                    j=1;
                                    while j<=jmax % loop 24
                                        jj=j+jump;
                                        if not(qw(jj)>0) % test 23
                                            hj=h*j;
                                            if j>jmax, qw(jj+1)=-1;end
                                            ehj=exp(hj);
                                            qw(jj)=q*ehj;
                                            vw(jj)= v*(hj+0.5-ehj*ehj*0.5);
                                        end % test 23
                                        pj=0;
                                        x=1-normcdf(gk-qw(jj))-pz;
                                        if x>0, pj=exp(w0+vw(jj)+r1*log(x)); end
                                        pk=pk+pj;
                                        if not(pj>pcutj) % test 24
                                            if jj>jmin || k>kmin % if conditions are met, make sure to exit loop 24
                                                j=jmax+1;
                                            end
                                        end
                                        j=j+1;
                                    end % loop 24
                                    h=-h;
                                    notDone22=h<0;
                                end % loop 22
                            end
                        end
                        prtrng=prtrng+pk;
                        if k>kmin && pk<=pcutk && pk1<=pcutk % if the conditions are met, the program will step out of loop 28 and loop 21
                            k=kmax+1; % make sure to exit loop 28
                            notDone21=false; % make sure to exit loop 21
                        else
                            pk2=pk1;
                            pk1=pk;
                            notDone21=gstep>0;
                        end
                    end % loop 21
                    k=k+1;
                end % loop 28
            end % if not(q<=0 || ifault==1)
        end % cdfTukey
        
    end % static methods
    
    methods
        
        function obj=StudRng(varargin)
            obj=obj@dContinuous('StudRng');
            obj.mindf = 4;
            obj.minr = 3;
            obj.CDFNearlyOne = 0.995;
            switch nargin
                case 0
                case 2
                    obj.NDistParms = 2;
                    obj.ParmTypes = 'ii';
                    obj.DefaultParmCodes = 'ii';
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('StudRng:Constructor', ...
                        'StudRng constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.df = newparmvalues(1);
            obj.r = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newdf  = ifelse(ParmCodes(1)=='f', obj.df, obj.df + 2);
            newr = ifelse(ParmCodes(2)=='f', obj.r, obj.r + 1);
            obj.ResetParms([newdf newr]);
        end
        
        function []=ReInit(obj)
            assert(obj.df>=obj.mindf,['StudRng df must be >= ' num2str(obj.mindf) '.']);
            assert(obj.r >=obj.minr, ['StudRng r must be >= '  num2str(obj.minr)  '.']);
            obj.Initialized = true;
            obj.LowerBound = 0.001;
            obj.UpperBound = 1000;
            % obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(obj.mindf,Parms(1))  NumTrans.GT2Real(obj.minr,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(obj.mindf,Reals(1)) NumTrans.Real2GT(obj.minr,Reals(2))];
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    thiscdf(iel) = obj.cdfTukey(X(iel),obj.df,obj.r);
                end
            end
        end
        
    end  % methods
    
end  % class StudRng

