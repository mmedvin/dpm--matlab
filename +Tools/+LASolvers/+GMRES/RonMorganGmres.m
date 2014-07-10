classdef RonMorganGmres < Tools.LASolvers.AbstractLA_Solver
%https://bearspace.baylor.edu/Ronald_Morgan/www/papers/gshift.pdf
    properties(Access = protected)
        A;
        %M;
        L;
        U;
        P;
        %ss;
        %ww;
        %p;
        
        Size;
        MaxSubSpace;
        SavekEVs;
        
        rtol;
        rtolev;
        cyclim;

        kperm;
        
        
        vk;
        hk;
        
        IsFirstTime=true;
        
        
        %debug
        %mvpl;
        mvp=0;
        gdr;
    end
    
    methods (Access = public)
        function obj = RonMorganGmres(A,params)
            obj.A = A;
            obj.Size = size(A,1);
            obj.MaxSubSpace = params.MaxSubSpace;
            obj.SavekEVs    = params.SavekEVs;
            
            obj.rtol = params.rtol;
            obj.rtolev = params.rtolev;
            obj.cyclim = params.cyclim;
            %obj.L = params.L;
            %obj.U = params.U;
            %obj.P = params.P;
            
            %
            obj.kperm = obj.SavekEVs;
            
            
        end
        
        function delete(obj)
%             figure(1)
%             clf reset
%             semilogy(obj.gdr(:,1));
%             
%             mvp1
             figure;
             clf reset
            semilogy(obj.gdr);
            title(num2str(obj.Size));
        end
        
        function x = FirstTimeSolve(obj,Rhs,x)
            m = obj.MaxSubSpace;
            k = obj.SavekEVs;
            
            
            %x = zeros(obj.Size,nrhs);
            %x = zeros(obj.Size,1);
            cycle = 1;
          %  mvp = 0;
            j = 1;
            rninit = norm(Rhs);
            rn = rninit;
            r = obj.minv(Rhs);
            vn = norm(r);
            v(:,1) = r/vn;
            %vn
            c(1,1) = vn;
            c(2:m+1,1) = zeros(m,1);
            
            h=zeros(m+1,m);
            %rna = zeros(obj.cyclim,m);
            %tha = zeros(obj.cyclim,m);
            %rhoa = zeros(obj.cyclim,m);
            rna = zeros(1,m);
            %tha = zeros(1,m);
            %rhoa = zeros(1,m);
            
            
            while ( (rn/rninit > obj.rtol) && (cycle <= obj.cyclim) )
                while ( j <= m )
                    wv = obj.op(v(:,j));
                    f = obj.minv(wv);
                    obj.mvp = obj.mvp + 1;
                    vnf = norm(f);
                    for i = 1:j
                        h(i,j) = v(:,i)'*f;
                        f = f - h(i,j) * v(:,i);
                    end  %for
                    vn = norm(f);
                    %--------------------------------------------------%
                    %--------reorthogonalization section-------%
                    if( vn < 0*vnf )   % A zero means reorthog is turned off.
                        % disp( 'did a reorthog')
                        for i = 1:j
                            dot_ = v(:,i)'*f;
                            f = f - dot_ * v(:,i);
                            h(i,j) = h(i,j) + dot_;
                        end  %for
                        vn = norm(f);
                    end  %if
                    %--------------------------------------------------%
                    h(j+1,j) = vn;
                    v(:,j+1) = f/h(j+1,j);
                    %subsection for res.norm at every iteration
                    d = h(1:j+1,1:j) \ c(1:j+1) ;
                    srv = c(1:j+1)-h(1:j+1,1:j)*d;
                    ojb.gdr(obj.mvp) = norm(srv);
                    %now subsection for res.norm if preconditioning used
                    %     xw = x(:,1) + v(:,1:j)*d;
                    %     wv = b(:,1) - op(n,xw);   %--may remove this and next line----%
                    %     gdr(mvp) = norm(wv);
                    %end subsection
                    j = j + 1;
                end  %while
                %hh = h(1:m,1:m);
                
                %----Set up and solve linear equations.-----%
                %  Should already be done:  d = h \ c;
                %
                %FOM section.
                %Next line is for FOM instead of GMRES:
                %  d = h(1:m,1:m) \ c(1:m);
                %End FOM section
                %
                %  srv = c-h*d;
                %cycle
                x = x + v(:,1:m)*d;
                %wv = Rhs(:,1) - obj.op(x(:,1)); %-may remove this and next line%
                %rnale(cycle) = norm(wv);
                %  r = minv(n,wv); etc.
                r = v(:,1:m+1)*srv;
                rn = norm(r);
                %residualnormforlinearequations = rn;
                %
                hh = h(1:m,1:m);
                em = zeros(m,1);
                em(m) = 1;
                ff = hh' \ em;
                %
                %FOM section
                %Comment next line for FOM instead of GMRES:
                hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
                %[g,dd] = eigs(hh,k,'sm');
                %[g,dd] = eig(full(hh),'nobalance');
                %ddd = diag(dd);
                %dabs = abs(ddd);
                %[thabs,ind] = sort(dabs);
                %th = ddd(ind);
                %ksmallestapproxeigenvalues = th(1:k);        % This prints the smallest k harmonic Ritz values.
                %gg = g(:,ind(1:k));
                
                [gg,dd] = eigs(hh,k,'sm');
                th=diag(dd);
                rho=zeros(k,1);
                for i=1:k
                    rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i);
                    tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
                    tvn = norm(tv);
                    %rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
                    %tha(cycle,i) = th(i);
                    %rhoa(cycle,i) = rho(i);
                    
                    rna(i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
                    %tha(i) = th(i);
                    %rhoa(i) = rho(i);

                    
                end  %for
                %
                greal = gg;
                i = 1;
                while( i <= k )
                    if( imag(th(i)) ~= 0 )
                        if( i ~= k )
                            greal(:,i+1) = imag(gg(:,i));
                            greal(:,i) = real(gg(:,i));
                            %disp('split complex vector')
                            i = i + 1;
                        else
                            k = k - 1;
                        end  %if
                    end  %if
                    i = i + 1;
                end  %while
                %
                greal(m+1,1:k) = zeros(1,k);
                greal(:,k+1) =   srv;
                [gon,rr] = qr(greal(:,1:k+1),0);
                hnew = gon'*h*gon(1:m,1:k);
                h(obj.kperm+1,:) = zeros(1,m);
                %
                %  Section for locking converged eigenvectors, right now it locks at
                %   rtolev*.1.  We don't really need to lock anymore, now that we are
                %   using 'nobalance' in the eig command.
                i = 1;
                while( i <= k )
                    if( rna(i) <= obj.rtolev*1.e-1 )
                        %if( rna(cycle,i) <= obj.rtolev*1.e-1 )
                        if( imag(th(i)) ~= 0 )
                            hnew(i+2:k+1,i) = zeros(k-i,1);
                            hnew(i+2:k+1,i+1) = zeros(k-i,1);
                            i = i + 1;
                        else
                            hnew(i+1:k+1,i) = zeros(k-i+1,1);
                        end  %if
                    end  %if
                    i = i + 1;
                end  %while
                %   end locking section
                %
                h(1:k+1,1:k) = hnew;
                c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
                c(k+2:m+1,1) = zeros(m-k,1);
                work = v*gon;
                v(:,1:k+1) = work;
                %
                %section for just reorthog. one vector, v_{k+1}
                for i = 1:k
                    dot_ = v(:,i)'*v(:,k+1) ;
                    v(:,k+1) = v(:,k+1) - dot_ * v(:,i);
                end  %for
                v(:,k+1) = v(:,k+1)/norm(v(:,k+1));
                % end section
                %
                j = k + 1;
                %kold = k;
                k = obj.kperm;
                cycle = cycle + 1;
            end  %while
            
            %obj.mvp1 = obj.mvp;
            
            k = obj.SavekEVs;% kold;
            obj.vk = v(:,1:k+1);
            obj.hk = h(1:k+1,1:k);
                        
        end
        
        function x = Solve(obj,Rhs,x)
            %
            % GPRO uses GMRES-DR for the first right-hand side followed by GMRES-Proj
            % for the other rhs's.  The 'd' at the end of 'gpro' means this is a copy
            % for distribution.  (But that does not mean it is efficiently programmed
            % or easy to understand!)
            %
            % subprogram 'op' is for the matrix multiply
            % subprogram 'minv' is for the preconditioner, but currently is the
            %     identity (it may or may not work for other preconditioners, sorry).
            %
            % n is the size of the matrix.
            % m is the max GMRES subspace.
            % k is the number of approximate eigenvectors saved at the restart.
            % rtol is the relative tolerance for solving the linear equations.
            % rtolev is sort of a tolerance for the e.vals that is used in locking them.
            % cyclim is the maximum number of GMRES-DR cycles.
            % nrhs is the number of right-hand sides
            %
            % The reorthogonalization section is currently turned off, but can be
            %   turned on by changing the zero on line 73 to 2.0 (or anything else > 1).
            % The value of 'm' used for GMRES in GMRES-Proj is currently set to m-k
            %   (for the m and k used for GMRES-DR).
            %
            % mvp is the number of matrix-vector products
            % gdr contains the residual norms at each iteration
            %
            %clear all;
            %format long e
            %global A
            %load sher4;
            %
            %n = 1104
            %Size = 2000; was n
            %m = 35;
            m = obj.MaxSubSpace;
            k = obj.SavekEVs;            
            
            if obj.IsFirstTime
                x = FirstTimeSolve(obj,Rhs,x);
                obj.IsFirstTime = false;
            else
                %  Section for solving other right-hand sides.
                %  This time the method is to alternate regular GMRES
                %    cycles with a projection over the already computed eigenvectors.
                %
                %  Next lines are in case there are changes in rtol and m and cyclim.
                %  As mentioned earlier, currently m is set to m-k.
                % rtol = 1.e-6
                % m = 15
               % m = m - obj.kperm;
                % cyclim = 50
                %
                
                %x = zeros(obj.Size,1);
                
                %h(1:m+1,1:m) = zeros(m+1,m);
                h = zeros(m+1,m);
                v = zeros(obj.Size,m);
                
                %for irhs = 2:nrhs
                %righthandside = irhs;
                
                cycle = 1;
                j = 1;
                rninit = norm(Rhs);
                rn = rninit;
                r = Rhs;
                while ( (rn/rninit > obj.rtol) && (cycle <= obj.cyclim) )
                    %
                    %  Projection section
                    %
                    iproj = 1;   %If iproj=1, then do the projection, otherwise it is turned off.
                    ifreq = 1;   %Frequency of projection over approx e.vects; currently projects between every GMRES cycle.
                    if( iproj == 1 )
                        if( cycle ==  1 + floor((cycle-1)/ifreq) * ifreq )
                            %if( cycle ==  floor(cycle/ifreq) * ifreq )
                            c = obj.vk'*r;     %This and next line for minres proj.
                            d = obj.hk  \ c ;
                            %     c = vk(:,1:k)'*r;  %This and next line for Galerkin proj.
                            %     d = hk(1:k,1:k) \ c;
                            %srv = c-hk*d;
                            x = x + obj.vk(:,1:k)*d;
                            wv = Rhs - obj.op(x); %-may remove this and next line%
                            obj.gdr(obj.mvp) = norm(wv);
                            r = obj.minv(wv);
                        end  %if involving ifreq
                    end  %if for iproj
                    %
                    %  end projection section
                    %
                    vn = norm(r);
                    v(:,1) = r/vn;
                    c(1,1) = vn;
                    c(2:m+1,1) = zeros(m,1);
                    %while ( (j <= m) ) %& (rn/rninit > rtol) )
                    for j=1:m
                        wv = obj.op(v(:,j));
                        f = obj.minv(wv);
                        obj.mvp = obj.mvp + 1;
                        vnf = norm(f);
                        for i = 1:j
                            h(i,j) = v(:,i)'*f;
                            f = f - h(i,j) * v(:,i);
                        end  %for
                        vn = norm(f);
                        %--------------------------------------------------%
                        %--------reorthogonalization section-------%
                        if( vn < 0*vnf )
                            % disp( 'did a reorthog')
                            for i = 1:j
                                dot = v(:,i)'*f;
                                f = f - dot * v(:,i);
                                h(i,j) = h(i,j) + dot;
                            end  %for
                            vn = norm(f);
                        end  %if
                        %--------------------------------------------------%
                        h(j+1,j) = vn;
                        v(:,j+1) = f/h(j+1,j);
                        %eyebar = eye(j+1,j);
                        %subsection for res.norm at every iteration
                        d = h(1:j+1,1:j) \ c(1:j+1) ;
                        srv = c(1:j+1)-h(1:j+1,1:j)*d;
                        %	srv
                        obj.gdr(obj.mvp) = norm(srv);
                        %rn = norm(srv);%gdr(mvp);
                        %now subsection for res.norm if preconditioning used
                        %     xw = x(:,irhs) + v(:,1:j)*d;
                        %     wv = b(:,irhs) - op(n,xw);   %--may remove this and next line----%
                        %     gdr(mvp) = norm(wv);
                        %end subsection
                        %j = j + 1;
                    end  %while
                    % if( rn/rninit > rtol )
                    %----Set up and solve linear equations.-----%
                    %
                    %FOM section.
                    %Next line is for FOM instead of GMRES:
                    %  d = h(1:m,1:m) \ c(1:m);
                    %End FOM section
                    %
                    %  srv = c-h*d;
                    %cycle
                    x = x + v(:,1:m)*d;
                    wv = Rhs - obj.op(x); %-may remove this and next line%
                    %rnale(cycle) = norm(wv);
                    %  r = minv(n,wv); etc.
                    % Note could change here :
                    %     r = v(:,1:m+1)*srv;
                    r = obj.minv(wv);
                    rn = norm(r);
                    %disp(rn)
                    j = 1;
                    cycle = cycle + 1;
                    % end  % if rtol...
                end  %while
            end
            
        end
        
        
         function x = SolveOld(obj,Rhs)            
%             %
%             % GPRO uses GMRES-DR for the first right-hand side followed by GMRES-Proj
%             % for the other rhs's.  The 'd' at the end of 'gpro' means this is a copy
%             % for distribution.  (But that does not mean it is efficiently programmed
%             % or easy to understand!)
%             %
%             % subprogram 'op' is for the matrix multiply
%             % subprogram 'minv' is for the preconditioner, but currently is the
%             %     identity (it may or may not work for other preconditioners, sorry).
%             %
%             % n is the size of the matrix.
%             % m is the max GMRES subspace.
%             % k is the number of approximate eigenvectors saved at the restart.
%             % rtol is the relative tolerance for solving the linear equations.
%             % rtolev is sort of a tolerance for the e.vals that is used in locking them.
%             % cyclim is the maximum number of GMRES-DR cycles.
%             % nrhs is the number of right-hand sides
%             %
%             % The reorthogonalization section is currently turned off, but can be
%             %   turned on by changing the zero on line 73 to 2.0 (or anything else > 1).
%             % The value of 'm' used for GMRES in GMRES-Proj is currently set to m-k
%             %   (for the m and k used for GMRES-DR).
%             %
%             % mvp is the number of matrix-vector products
%             % gdr contains the residual norms at each iteration
%             %
%             %clear all;
%             %format long e
%             %global A
%             %load sher4;
%             %
%             %n = 1104
%             %Size = 2000; was n
%             %m = 35; 
%             m = obj.MaxSubSpace;
%             %k = 15;
%             k = obj.SavekEVs;
%             %rtol = 1.e-8;
%             %rtolev = 1.e-16;
%             %cyclim = 30;
%             nrhs = size(Rhs,2);
%             %
%             %kperm = obj.SavekEVs;
%             %---- Matrix can be set up here or can be set up in subprogram 'op'. ----%
%             %rand('seed',0)
%             %A = spdiags(sort(rand(Size,1)),0,Size,Size) + spdiags(.0002*ones(Size,1),1,Size,Size);
%             %-------------- RHS vectors.--------------------------%
%             %randn('seed',0)
%             %Rhs = randn(Size,nrhs);
%             %-------------------------------------------------------%
%             x = zeros(obj.Size,nrhs);
%             cycle = 1;
%             mvp = 0;
%             j = 1;
%             rninit = norm(Rhs(:,1));
%             rn = rninit;
%             r = obj.minv(Rhs(:,1));
%             vn = norm(r);
%             v(:,1) = r/vn;
%             %vn
%             c(1,1) = vn;
%             c(2:m+1,1) = zeros(m,1);
%             while ( (rn/rninit > obj.rtol) && (cycle <= obj.cyclim) )
%                 while ( j <= m )
%                     wv = obj.op(v(:,j));
%                     f = obj.minv(wv);
%                     mvp = mvp + 1;
%                     vnf = norm(f);
%                     for i = 1:j
%                         h(i,j) = v(:,i)'*f;
%                         f = f - h(i,j) * v(:,i);
%                     end  %for
%                     vn = norm(f);
%                     %--------------------------------------------------%
%                     %--------reorthogonalization section-------%
%                     if( vn < 0*vnf )   % A zero means reorthog is turned off.
%                         % disp( 'did a reorthog')
%                         for i = 1:j
%                             dot = v(:,i)'*f;
%                             f = f - dot * v(:,i);
%                             h(i,j) = h(i,j) + dot;
%                         end  %for
%                         vn = norm(f);
%                     end  %if
%                     %--------------------------------------------------%
%                     h(j+1,j) = vn;
%                     v(:,j+1) = f/h(j+1,j);
%                     %subsection for res.norm at every iteration
%                     d = h(1:j+1,1:j) \ c(1:j+1) ;
%                     srv = c(1:j+1)-h(1:j+1,1:j)*d;
%                     gdr(mvp,1) = norm(srv);
%                     %now subsection for res.norm if preconditioning used
%                     %     xw = x(:,1) + v(:,1:j)*d;
%                     %     wv = b(:,1) - op(n,xw);   %--may remove this and next line----%
%                     %     gdr(mvp) = norm(wv);
%                     %end subsection
%                     j = j + 1;
%                 end  %while
%                 %hh = h(1:m,1:m);
%                 
%                 %----Set up and solve linear equations.-----%
%                 %  Should already be done:  d = h \ c;
%                 %
%                 %FOM section.
%                 %Next line is for FOM instead of GMRES:
%                 %  d = h(1:m,1:m) \ c(1:m);
%                 %End FOM section
%                 %
%                 %  srv = c-h*d;
%                 %cycle
%                 x(:,1) = x(:,1) + v(:,1:m)*d;
%                 wv = Rhs(:,1) - obj.op(x(:,1)); %-may remove this and next line%
%                 rnale(cycle) = norm(wv);
%                 %  r = minv(n,wv); etc.
%                 r = v(:,1:m+1)*srv;
%                 rn = norm(r);
%                 %residualnormforlinearequations = rn;
%                 %
%                 hh = h(1:m,1:m);
%                 em = zeros(m,1);
%                 em(m) = 1;
%                 ff = hh' \ em;
%                 %
%                 %FOM section
%                 %Comment next line for FOM instead of GMRES:
%                 hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
%                 %[g,dd] = eigs(hh,k);
%                 [g,dd] = eig(full(hh),'nobalance');
%                 ddd = diag(dd);
%                 dabs = abs(ddd);
%                 [thabs,ind] = sort(dabs);
%                 th = ddd(ind);
%                 %ksmallestapproxeigenvalues = th(1:k);        % This prints the smallest k harmonic Ritz values.
%                 gg = g(:,ind(1:k));
%                 for i=1:k
%                     rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i);
%                     tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
%                     tvn = norm(tv);
%                     rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
%                     tha(cycle,i) = th(i);
%                     rhoa(cycle,i) = rho(i);
%                 end  %for
%                 %
%                 greal = gg;
%                 i = 1;
%                 while( i <= k )
%                     if( imag(th(i)) ~= 0 )
%                         if( i ~= k )
%                             greal(:,i+1) = imag(gg(:,i));
%                             greal(:,i) = real(gg(:,i));
%                             %disp('split complex vector')
%                             i = i + 1;
%                         else
%                             k = k - 1;
%                         end  %if
%                     end  %if
%                     i = i + 1;
%                 end  %while
%                 %
%                 greal(m+1,1:k) = zeros(1,k);
%                 greal(:,k+1) =   srv;
%                 [gon,rr] = qr(greal(:,1:k+1),0);
%                 hnew = gon'*h*gon(1:m,1:k);
%                 h(obj.kperm+1,:) = zeros(1,m);
%                 %
%                 %  Section for locking converged eigenvectors, right now it locks at
%                 %   rtolev*.1.  We don't really need to lock anymore, now that we are
%                 %   using 'nobalance' in the eig command.
%                 i = 1;
%                 while( i <= k )
%                     if( rna(cycle,i) <= obj.rtolev*1.e-1 )
%                         if( imag(th(i)) ~= 0 )
%                             hnew(i+2:k+1,i) = zeros(k-i,1);
%                             hnew(i+2:k+1,i+1) = zeros(k-i,1);
%                             i = i + 1;
%                         else
%                             hnew(i+1:k+1,i) = zeros(k-i+1,1);
%                         end  %if
%                     end  %if
%                     i = i + 1;
%                 end  %while
%                 %   end locking section
%                 %
%                 h(1:k+1,1:k) = hnew;
%                 c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
%                 c(k+2:m+1,1) = zeros(m-k,1);
%                 work = v*gon;
%                 v(:,1:k+1) = work;
%                 %
%                 %section for just reorthog. one vector, v_{k+1}
%                 for i = 1:k
%                     dot = v(:,i)'*v(:,k+1) ;
%                     v(:,k+1) = v(:,k+1) - dot * v(:,i);
%                 end  %for
%                 v(:,k+1) = v(:,k+1)/norm(v(:,k+1));
%                 % end section
%                 %
%                 j = k + 1;
%                 kold = k;
%                 k = obj.kperm;
%                 cycle = cycle + 1;
%             end  %while
%             
%            % mvp1 = mvp
%             %
%             k = kold;
%             vk = v(:,1:k+1);
%             hk = h(1:k+1,1:k);
%             %
%             %
%            % figure(1)
%            % clf reset
%            % semilogy(gdr(:,1));
%             %
%             %  Section for solving other right-hand sides.
%             %  This time the method is to alternate regular GMRES
%             %    cycles with a projection over the already computed eigenvectors.
%             %
%             %  Next lines are in case there are changes in rtol and m and cyclim.
%             %  As mentioned earlier, currently m is set to m-k.
%             % rtol = 1.e-6
%             % m = 15
%             m = m - obj.kperm;
%             % cyclim = 50
%             %
%             h(1:m+1,1:m) = zeros(m+1,m);
%             for irhs = 2:nrhs
%                 %righthandside = irhs;
%                 cycle = 1;
%                 j = 1;
%                 rninit = norm(Rhs(:,irhs));
%                 rn = rninit;
%                 r = Rhs(:,irhs);
%                 while ( (rn/rninit > obj.rtol) && (cycle <= obj.cyclim) )
%                     %
%                     %  Projection section
%                     %
%                     iproj = 1;   %If iproj=1, then do the projection, otherwise it is turned off.
%                     ifreq = 1;   %Frequency of projection over approx e.vects; currently projects between every GMRES cycle.
%                     if( iproj == 1 )
%                         if( cycle ==  1 + floor((cycle-1)/ifreq) * ifreq )
%                             %if( cycle ==  floor(cycle/ifreq) * ifreq )
%                             c = vk'*r;     %This and next line for minres proj.
%                             d = hk  \ c ;
%                             %     c = vk(:,1:k)'*r;  %This and next line for Galerkin proj.
%                             %     d = hk(1:k,1:k) \ c;
%                             %srv = c-hk*d;
%                             x(:,irhs) = x(:,irhs) + vk(:,1:k)*d;
%                             wv = Rhs(:,irhs) - obj.op(x(:,irhs)); %-may remove this and next line%
%                             gdr(mvp) = norm(wv);
%                             r = obj.minv(wv);
%                         end  %if involving ifreq
%                     end  %if for iproj
%                     %
%                     %  end projection section
%                     %
%                     vn = norm(r);
%                     v(:,1) = r/vn;
%                     c(1,1) = vn;
%                     c(2:m+1,1) = zeros(m,1);
%                     while ( (j <= m) ) %& (rn/rninit > rtol) )
%                         wv = obj.op(v(:,j));
%                         f = obj.minv(wv);
%                         mvp = mvp + 1;
%                         vnf = norm(f);
%                         for i = 1:j
%                             h(i,j) = v(:,i)'*f;
%                             f = f - h(i,j) * v(:,i);
%                         end  %for
%                         vn = norm(f);
%                         %--------------------------------------------------%
%                         %--------reorthogonalization section-------%
%                         if( vn < 0*vnf )
%                             % disp( 'did a reorthog')
%                             for i = 1:j
%                                 dot = v(:,i)'*f;
%                                 f = f - dot * v(:,i);
%                                 h(i,j) = h(i,j) + dot;
%                             end  %for
%                             vn = norm(f);
%                         end  %if
%                         %--------------------------------------------------%
%                         h(j+1,j) = vn;
%                         v(:,j+1) = f/h(j+1,j);
%                         eyebar = eye(j+1,j);
%                         %subsection for res.norm at every iteration
%                         d = h(1:j+1,1:j) \ c(1:j+1) ;
%                         srv = c(1:j+1)-h(1:j+1,1:j)*d;
%                         %	srv
%                         gdr(mvp) = norm(srv);
%                         rn = gdr(mvp);
%                         %now subsection for res.norm if preconditioning used
%                         %     xw = x(:,irhs) + v(:,1:j)*d;
%                         %     wv = b(:,irhs) - op(n,xw);   %--may remove this and next line----%
%                         %     gdr(mvp) = norm(wv);
%                         %end subsection
%                         j = j + 1;
%                     end  %while
%                     % if( rn/rninit > rtol )
%                     %----Set up and solve linear equations.-----%
%                     %
%                     %FOM section.
%                     %Next line is for FOM instead of GMRES:
%                     %  d = h(1:m,1:m) \ c(1:m);
%                     %End FOM section
%                     %
%                     %  srv = c-h*d;
%                     %cycle
%                     x(:,irhs) = x(:,irhs) + v(:,1:m)*d;
%                     wv = Rhs(:,irhs) - obj.op(x(:,irhs)); %-may remove this and next line%
%                     rnale(cycle) = norm(wv);
%                     %  r = minv(n,wv); etc.
%                     % Note could change here :
%                     %     r = v(:,1:m+1)*srv;
%                     r = wv;
%                     rn = norm(r);
%                     %disp(rn)
%                     j = 1;
%                     cycle = cycle + 1;
%                     % end  % if rtol...
%                 end  %while
%                 %mvp
%             end  %for each rhs
%             %mvp1
%             %figure(2)
%             %clf reset
%             %semilogy(gdr(:,1));
%             
%             
         end
        
    end
    
    methods (Access = protected)
        
        function [y] = op(obj,x)
            %global A
            %global ss
            %global ww
            %global p
            y = obj.A*x;            
            % y = y + ss.*x;
            % y = y + (ww'*x)*ss;
        end
        
        function [y] = minv(obj,x)
            %global D
            % Diagonal preconditioning:
            %y = D.\x;
            %Preconditioner M
            %global M
            %y = obj.U\(obj.L\(obj.P*x));
            % % Next line is for no preconditioning.
            y = x;
        end
        
    end
    
end

