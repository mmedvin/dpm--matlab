%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function demo1_run

demo_globals

az = 130; el = 40;

nx = 23;  h = 1/(nx+1);

[A, N] = sp_laplace(nx,nx);

rand('seed', 1234567);      xt = 2 * rand(N,1) - 1;

b = A*xt;
x     =  zeros(N,1);
error =  xt - x;
ERROR = reshape(error, nx, nx);

subplot(1,1,1)
cla
subplot(2,2,1)
hold off
surf(ERROR)
hold on
title(' Initial error, Physical Space')
hold off
ERROR_COEF = sint2(ERROR);
subplot(2,2,2)
hold off
surf(ERROR_COEF); view([az,el])
hold on
title(' Initial error, Fourier Space')
hold off
pause(4)

smooth = DEMO_VAR1;

if smooth == 1
   % weighted Jacobi
   D = 0.95 * diag(diag(A));
elseif smooth == 2
   % Gauss-Seidel
   L = tril(A);
elseif smooth == 3
   % red/black Gauss-Seidel
   red = [1:2:N]; black = [2:2:N];
   NR=length(red); NB = length(black);
   p = [red, black];
   ip = zeros(N,1);
   for i = 1:N, ip(p(i))=i; end
   R = A(red,red); B = A(black,black); C = A(red,black);
   b_p = b(p); br = b_p(1:NR); bb = b_p(NR+1:N);
   xr = zeros(NR,1); xb = zeros(NB,1);
end
 
j = 0;
while j < 6
   j = j+1;
   if smooth == 1
      x = x + D \ (b - A*x);
   elseif smooth == 2
      x = x + L \(b - A*x);
   elseif smooth == 3
      xr = R \ (br - C  * xb);
      xb = B \ (bb - C' * xr);
      x_p = [xr; xb];
      x = x_p(ip);
   end
 
   error =  xt - x;
   ERROR = reshape(error, nx, nx);
   subplot(2,2,1)
   hold off
   surf(ERROR)
   hold on
   if smooth == 1
   title(['Physical Space, ',num2str(j),' iter(s), weighted Jacobi.'])
   elseif smooth == 2
   title(['Physical Space, ',num2str(j),' iter(s), Gauss-Seidel.'])
   elseif smooth == 3
   title(['Physical Space, ',num2str(j),' iter(s), Red-Black Gauss-Seidel.'])
   end

   hold off
   ERROR_COEF = sint2(ERROR);
   subplot(2,2,2)
   hold off
   surf(ERROR_COEF); view([az,el])
   hold on
   if smooth == 1
   title(['Fourier Space, ',num2str(j),' iter(s), weighted Jacobi.'])
   elseif smooth == 2
   title(['Fourier Space, ',num2str(j),' iter(s), Gauss-Seidel.'])
   elseif smooth == 3
   title(['Fourier Space, ',num2str(j),' iter(s), R/B Gauss-Seidel.'])
   end

   hold off

   j = j+1;
   if smooth == 1
      x = x + D \ (b - A*x);
   elseif smooth == 2
      x = x + L \(b - A*x);
   elseif smooth == 3
      xr = R \ (br - C  * xb);
      xb = B \ (bb - C' * xr);
      x_p = [xr; xb];
      x = x_p(ip);
   end
 
   error =  xt - x;
   ERROR = reshape(error, nx, nx);
   subplot(2,2,3)
   hold off
   surf(ERROR)
   hold on
   if smooth == 1
   title(['Physical Space, ',num2str(j),' iter(s), weighted Jacobi.'])
   elseif smooth == 2
   title(['Physical Space, ',num2str(j),' iter(s), Gauss-Seidel.'])
   elseif smooth == 3
   title(['Physical Space, ',num2str(j),' iter(s), Red-Black Gauss-Seidel.'])
   end

   hold off
   ERROR_COEF = sint2(ERROR);
   subplot(2,2,4)
   hold off
   surf(ERROR_COEF); view([az,el])
   hold on
   if smooth == 1
   title(['Fourier Space, ',num2str(j),' iter(s), weighted Jacobi.'])
   elseif smooth == 2
   title(['Fourier Space, ',num2str(j),' iter(s), Gauss-Seidel.'])
   elseif smooth == 3
   title(['Fourier Space, ',num2str(j),' iter(s), Red-Black Gauss-Seidel.'])
   end

   hold off

   pause(4)
end

hold off

% Reset axes
% close(fig1); close; MGLab
