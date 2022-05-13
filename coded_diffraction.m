clear;
close all;
clc;

tic

RGB_img = im2double(imread('venere.jpg'));

n1 = size(RGB_img, 1);
n2 = size(RGB_img, 2);
d = n1*n2;
Lgrid = 1 : 1 : 20; % number of patterns
Tgrid = 10000; % number of iterations of power method
M = 40;
niter = length(Lgrid);
Delta = 10;
eps = 10^(-7);
deltatrg = 1.001;

scal = zeros(niter, max(Tgrid), 3);

x0hat = zeros(d, 3, niter);

for iter = 1 : niter
    L = Lgrid(iter);
    alpha = 100*d;
    T = Tgrid;
    
    fprintf('L=%d\n', L);


for i = 1 : 3
   
   fprintf('Image %d\n', i);
   
   x = RGB_img(:, :, i);
   x = x/sqrt(sum(sum(x.^2)));
   xlin = reshape(x, 1, d);
   
   cdp = zeros(n1, n2, L);
   
   for j1 = 1 : L
       
       u = rand(1, d);
       b1 = ones(1, d);
       b1(u>1/4) = -1;
       b1(u>1/2) = -1i;
       b1(u>3/4) = 1i;
          
       cdp(:, :, j1) = reshape(b1, n1, n2);
   end
   
   % generation of observations
   ymat = zeros(n1, n2, L);
   
   for j1 = 1 : L
       xmod = x .* conj(cdp(:, :, j1));
       ymat(:, :, j1) = 1-sqrt(deltatrg)./(abs(fft2(xmod)).^2-1+sqrt(deltatrg));
       ymat(ymat<-M)=-M;

   end
   
   % power method
   
   % initialization
   x0norm = 1/sqrt(2) * (randn(n1, n2) + 1i * randn(n1, n2));
   x0 = x0norm/sqrt(sum(sum(abs(x0norm).^2)));
   
   x0old = zeros(Delta, d);
   
   for jit = 1 : T
       
       mul1mat = zeros(n1, n2, L);

       for j1 = 1 : L
           x0mod = x0 .* conj(cdp(:, :, j1));
           mul1mat(:, :, j1) = fft2(x0mod);
       end
              
       mul2 = mul1mat .* ymat;
       
       mul3 = zeros(n1, n2, L);
       
       for j1 = 1 : L
           mul3(:, :, j1) = d*ifft2(mul2(:, :, j1)).* cdp(:, :, j1);
       end
       
       x0norm = sum(mul3, 3) + alpha * x0;
       x0 = x0norm/sqrt(sum(sum(abs(x0norm).^2))); 
       x0lin = reshape(x0, 1, d);
       
       fprintf('Iteration %d: %f', jit, abs(xlin * x0lin'));
       scal(iter, jit, i) = abs(xlin * x0lin');
       save data_2Dfft_optnew_venere.mat scal Lgrid;
       
       if jit <= Delta
           x0old(jit, :) = x0lin;
           fprintf('\n');
       else
           fprintf(', scalar with previous: %.8f\n', abs(x0old(1, :) * x0lin'));
           if abs(x0old(1, :) * x0lin') > 1- eps
               break;
           else
               for jjj = 1 : 9
                   x0old(jjj, :) = x0old(jjj+1, :);
               end
               x0old(10, :) = x0lin;
           end
       end                  

   end
   
   x0hat(:, i, iter) = x0lin;
   save data_2Dfft_optnew_venere.mat scal Lgrid x0hat;

   fprintf('\n');
   
end
   fprintf('\n');

end

toc