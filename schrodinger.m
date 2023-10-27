close all;

% Initialize variables
N = 256; L = 20; [D,x] = cheb(N);
x = x*L; D2 = D^2; D2 = D2(2:N,2:N);
[~,w] = clencurt(N); w = w*L;

% Initial condition
psi0 = exp(-0.1*x.^2); psi0(1)=0; psi0(N+1)=0;
psi0 = psi0/sqrt(w*psi0);
% psi0 = psi0/sqrt(nonUniformTrap(x,psi0.*psi0));

% Construct Hamiltonian operator
V = zeros(N+1); V = V(2:N); % Arbitrary potential
H = -0.5*D2 + diag(V);

% Compute eigenfunctions and eigenvalues, sort
[P,DD] = eig(H);
[E,ind] = sort(diag(DD));
P = P(:,ind);
%%
% Columns of P are eigenfunctions and E is a vector of energy eigenvalues.

% Normalize the eigenfunctions (each squared integrates to 1)
for i=1:N-1
    int_val = w * ([0;P(:,i);0].*[0;conj(P(:,i));0]);
    % int_val = nonUniformTrap(x,[0;P(:,i);0].*[0;conj(P(:,i));0]);
    % int_val = trapz(x,-[0;P(:,i);0]);
    P(:,i) = P(:,i)/sqrt(int_val);
end

%%
% Compute c_j coefficients
c = zeros(N-1,1);
for i=1:N-1
    c(i) = w * ([0;psi0(2:N);0] .* [0;conj(P(:,i));0]);
%     c(i) = dot([0;psi0(2:N);0], [0;conj(P(:,i));0]);
    % c(i)=nonUniformTrap(x,([0;psi0(2:N);0] .* [0;conj(P(:,i));0]));
end

% Construct solution
dt = 0.005; tf = 10; tvec = 0:dt:tf;
tsteps = tf/dt;
soln = zeros(N+1,tsteps+1);
soln(:,1) = psi0;
for i=1:tsteps
    psi = sum(P*(c.*exp(1i*E*tvec(i))), 2);
    soln(:,i+1) = [0;psi;0];
end
soln = soln';

%%

figure
hold on
for i=1:4
    plot(x,[0;P(:,i);0],LineWidth=2);
end
legend('$n=1$','$n=2$','$n=3$','$n=4$',Interpreter='latex')

% %%
% 
% figure
% hold on
% xx = linspace(-1,1,1000);
% uu = polyval(polyfit(x,psi0,N),xx);
% plot(xx,uu,'LineWidth',2);
% scatter(x,psi0);
% 
%%

figure
for i=1:tsteps+1
    plot(x, abs(soln(i,:)).^2, LineWidth=2);
    ylim([-.1,.1]);
    title(['$t=',num2str(tvec(i)),'$'],'Interpreter','latex')
    drawnow;
    pause(0.005);
end