k1=1;k2=2;
for Kind=3%1:2
    errpre=0;%
    
    for N = 1:9
        
        %create grid
        n=2^(N+1)+1;
        x=linspace(-1,1,n+2);
        y=linspace(-1,1,n+2);
        h= abs(x(1)-x(2));
        [X,Y]=meshgrid(x(2:end-1),y(2:end-1));
        
        %create descrete operator
        I = speye(n,n);
        E = sparse(2:n,1:n-1,1,n,n);
        D = E+E'-2*I;
        A = (kron(D,I)+kron(I,D))./(h^2);
        
        %create RHS and the exact solution
        if Kind==1
            f=2* ((X.^2 - 1) + (Y.^2 - 1));
            ex=(X.^2 - 1) .*(Y.^2 - 1);
        elseif Kind==2
            ex = sin(pi*X).*sin(pi*Y);
            f=-ex*2*pi^2;
        elseif Kind==3
            g=(k2-k1)./(1+exp(-10000*Y)) +k1;
            f = 2* ((X.^2 - 1) + (Y.^2 - 1))./g;
            %f = -2*pi^2* (sin(pi*X).*sin(pi*Y))./g;
            
        end
        %solve the problem
        ap = (A\f(:));
        
        
        if Kind==3
            if N>1
                tmp=reshape(ap,n,n);
                tmp=tmp(1:2:end,1:2:end);
                err  = norm(tmp(:)-v(:),inf)/norm(tmp(:),inf);

                %print the error and the rate
                fprintf('Kind=%d\t n=%d\t err=%-10.8d\t rate=%-4.2f\n',Kind,n,err,log2(errpre/err));
                
                %keep the old value for the rate calculation at next step
                errpre=err;

            end
            v=ap;
        else
            %compute the error
            err  = norm(ex(:)-ap(:),inf)/norm(ap(:),inf);
            
            %print the error and the rate
            fprintf('Kind=%d\t n=%d\t err=%-10.8d\t rate=%-4.2f\n',Kind,n,err,log2(errpre/err));
            
            %keep the old value for the rate calculation at next step
            errpre=err;
        end
       
        
    end
end
