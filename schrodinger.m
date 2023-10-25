% Initialize variables
N = 16; L = 1; [D,x] = cheb(N);
x = x*L; D2 = D^2; D2 = D2(2:N,2:N);
psi0 = exp(-x.^2);

% Construct Hamiltonian operator
V = zeros(N+1,1); V = V(2:N);
H = -0.5*D2 + diag(V);

% Compute eigenfunctions and eigenvalues
[P,DD] = eig(H); E = diag(DD);

% Compute c_j coefficients

% Construct solution

figure
hold on
for i=1:8
    plot(x,[0;P(:,i);0]);
end