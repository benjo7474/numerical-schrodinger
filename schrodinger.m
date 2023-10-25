% Initialize variables
N = 64; L = 20; [D,x] = cheb(N);
x = x*L; D2 = D^2; D2 = D2(2:N,2:N);
psi0 = exp(-0.1*x.^2); psi0 = psi0(2:N);

% Construct Hamiltonian operator
V = x.^2; V = V(2:N); % Arbitrary potential
H = -0.5*D2 + diag(V);

% Compute eigenfunctions and eigenvalues, sort
[P,DD] = eig(H);
[E,ind] = sort(diag(DD));
P = P(:,ind);
% Columns of P are eigenfunctions and E is a vector of energy eigenvalues.

% Normalize the eigenfunctions (each integrates to 1)
for i=1:N-1
    int_val = trapz(x, [0;P(:,i).*conj(P(:,i));0]);
    P(:,i) = P(:,i)/int_val;
end

% Compute c_j coefficients
c = zeros(N-1,1); % represents the probabilities of finding particle energy
for i=1:N-1
    c(i) = trapz(x, [0;psi0;0].*conj([0;P(:,i);0])); % Need to justify this
end

% Normalize the c_j coefficients
c = c/sum(c);

% Construct solution
dt = 0.01; tf = 1; tvec = 0:dt:tf;
tsteps = tf/dt;
soln = zeros(N+1,tsteps+1);

% Normalize the solution


%%

figure
hold on
for i=1:4
    plot(x,[0;P(:,i);0],LineWidth=2);
end
legend('$n=1$','$n=2$','$n=3$','$n=4$',Interpreter='latex')

%% 

figure
for i=1:tsteps+1
    plot(x,abs(soln(i,:)).^2,LineWidth=2)
    title(['$t=',num2str(tvec(i)),'$'],'Interpreter','latex')
    ylim([-2,2]);
    drawnow;
    pause(0.02)
end