% Initialize variables
N = 32; L = 10; [D,x] = cheb(N);
x = x*L; D2 = D^2; D2 = D2(2:N,2:N);
psi0 = exp(-x.^2); psi0 = psi0(2:N);

% Construct Hamiltonian operator
V = x.^2; V = V(2:N); % Arbitrary potential
H = -0.5*D2 + diag(V);

% Compute eigenfunctions and eigenvalues, sort
[P,DD] = eig(H);
[E,ind] = sort(diag(DD));
P = P(:,ind);

% Compute c_j coefficients and construct solution
c = zeros(N-1,1); % represents the probabilities of finding particle energy
dt = 0.01; tf = 1; tvec = 0:dt:tf;
tsteps = tf/dt;
soln = zeros(N+1,tsteps+1);
for i=1:N-1
    c(i) = sum(psi0.*P(:,i))/L; % Trapezoid rule
    soln = soln + c(i)*[0;P(:,i);0]*exp(1i*E(i)*tvec);
end
soln = soln';

% Normalize each row
for i=1:tsteps+1
    soln(i,:) = soln(i,:) / (abs(soln(i,:)).^2);
end

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
    ylim([-200,1000]);
    drawnow;
    pause(0.05)
end