clear;
close all;
clc;

tic

delta = 0.4:0.4:4;
d = 4096;
n = floor(d * delta);
numT = 40;
deltatrg = 1.001/2;

for i = 1 : length(n)
    fprintf('n=%d\n', n(i));
    for ii = 1 : numT
    
    
        fprintf('sample=%d\n', ii);
        x0norm = randn(d, 1);
        x0 = x0norm/sqrt(sum(x0norm.^2));
        A = randn(d, n(i));
        z = zeros(1, n(i));
        psi = zeros(1, n(i));

        for j = 1 : n(i)
            z(j) = sum(A(:, j) .* x0);
            psi(j) = 1-sqrt(2*deltatrg)/(z(j)^2-1+sqrt(2*deltatrg));                                    
        end
    
        M = 1/n(i) * A * diag(psi)* A';
    
        [Vc, Dc] = eig(M);
        Vr = real(Vc);
        Dr = real(Dc);
        
        [valeig, indeig]=sort(diag(Dr), 'descend');
        v0 = Vr(:, indeig(1));
    
        scal(ii, i) = (sum(v0.* x0))^2;
        save data_opt_2delta1dot001_d4096.mat n d delta scal;
    end
    fprintf('\n');
    
end

save data_opt_2delta1dot001_d4096.mat n d delta scal;

toc