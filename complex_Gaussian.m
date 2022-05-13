clear;
close all;
clc;

tic

delta = 0.8:0.8:8;
d = 4096;
n = floor(d * delta);
numT = 40;
deltatrg = 1.001;


for i1 = 1 : length(n)
    fprintf('n=%d\n', n(i1));
    for ii = 1 : numT
    
        fprintf('sample=%d\n', ii);
        x0norm = 1/sqrt(2) * (randn(d, 1) + 1i * randn(d, 1));
        x0 = x0norm/sqrt(sum(abs(x0norm).^2));
        A = 1/sqrt(2) * (randn(d, n(i1)) + 1i * randn(d, n(i1)));
        z = zeros(1, n(i1));
        psi = zeros(1, n(i1));

        for j = 1 : n(i1)
            z(j) = abs(sum(A(:, j) .* conj(x0)));
            psi(j) = 1-sqrt(deltatrg)/(z(j)^2-1+sqrt(deltatrg));                        
        end
    
        M = 1/n(i1) * A * diag(psi)* A';

        [Vc, Dc] = eig(M);
        
        [valeig, indeig]=sort(diag(real(Dc)), 'descend');

        v0 = Vc(:, indeig(1));    
    
        scal(ii, i1) = (abs(sum(v0.* conj(x0))))^2;
        save data_opt_delta1dot001_d4096.mat n d delta scal;
    end
    fprintf('\n');
    
end

save data_opt_delta1dot001_d4096.mat n d delta scal;

toc