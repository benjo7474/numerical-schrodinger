%% Solve system

clear; close all;

% Initialize variables
N = 256; L = 10; [D,x] = cheb(N);
D = D/L; x = x*L; D2 = D^2; D2 = D2(2:N,2:N); 
[~,w] = clencurt(N); w = w*L;

% Specify initial condition and potential function
% psi0 = exp(-(x+5).^2/2).*exp(-1i*1.*x); % IC
psi0 = exp(-(x).^2/2);
% V = 30*exp(-0.5*(x-10).^2) + 30*exp(-0.5*(x+10).^2); % V
% V = 0.5*x.^2;
V = zeros(N+1);
% V = 0.1*pot(x);

% Make sure IC is normalized
psi0 = psi0/sqrt(w*(abs(psi0).^2));
% psi0 = psi0/sqrt(nonUniformTrap(x,psi0.*psi0));

% Construct Hamiltonian operator
V = V(2:N);
H = -0.5*D2 + diag(V);

% Compute eigenfunctions and eigenvalues, sort by smallest eigenvalues
[P,DD] = eig(H);
[E,ind] = sort(diag(DD));
P = P(:,ind);
% Columns of P are eigenfunctions and E is a vector of energy eigenvalues.

% Normalize the eigenfunctions (each squared integrates to 1)
for i=1:N-1
    int_val = w * ([0;P(:,i);0].*[0;conj(P(:,i));0]);
    % int_val = nonUniformTrap(x,[0;P(:,i);0].*[0;conj(P(:,i));0]);
    % int_val = trapz(x,-[0;P(:,i);0]);
    P(:,i) = P(:,i)/sqrt(int_val);
end

% Compute c_j coefficients
c = zeros(N-1,1);
for i=1:N-1
    c(i) = w * ([0;psi0(2:N);0] .* [0;conj(P(:,i));0]);
    % c(i) = dot([0;psi0(2:N);0], [0;conj(P(:,i));0]);
    % c(i)=nonUniformTrap(x,([0;psi0(2:N);0] .* [0;conj(P(:,i));0]));
end

% Construct solution
dt = 0.05; tf = 100; tvec = 0:dt:tf;
tsteps = tf/dt;
soln = zeros(N+1,tsteps+1);
soln(:,1) = psi0;
for i=1:tsteps
    psi = sum(P*(c.*exp(1i*E*tvec(i))), 2);
    soln(:,i+1) = [0;psi;0];
end
soln = soln';

%% Defaults
set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',5);
set(0,'defaultAxesFontSize',35)

%% Plot first couple eigenfunctions (numerical)

figure
hold on

plot(x,[0;P(:,1);0],Color='#f2ce18');
plot(x,[0;P(:,2);0],Color='#38cfc2');
plot(x,[0;P(:,3);0],Color='#f218e0');
plot(x,[0;P(:,4);0],Color='r');

legend('$n=1$','$n=2$','$n=3$','$n=4$')
xlabel('$x$')
ylabel('$\psi_n(x)$')
title('Eigenfunctions for Infinite Square Well')

%% Plot initial condition
% 
% figure
% hold on
% xx = linspace(-1,1,1000);
% uu = polyval(polyfit(x,psi0,N),xx);
% plot(xx,uu,'LineWidth',2);
% scatter(x,psi0);

%% Plot potential
figure
plot(x(2:N), V, '-r', 'LineWidth', 2);
title('Potential', 'FontSize', 16)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$V(x)$', 'Interpreter', 'latex', 'FontSize', 16);

%% Plot eigenvalues
figure
n = 1:10;
xx = linspace(0,10.5,1000);
uu = pi^2/(8*L^2)*xx.^2;
hold on
plot(xx,uu,'Color','#38cfc2');
scatter(n, E(n), 200, 'k', 'filled');
legend('$n^2\pi^2/(8L^2)$','$E_n$ (numerical)','Location','northwest','FontSize',40)
title('Eigenvalues for Infinite Square Well')

%% Animiation of solution

figure
ymax = max(abs(soln).^2, [], 'all');
for i=1:tsteps+1
    plot(x, abs(soln(i,:)).^2, LineWidth=2);
    hold on
    plot(x(2:N), 0.01*V, '-r')
    hold off
    ylim([-0.1*ymax, 1.1*ymax]);
    xlabel('$x$')
    ylabel('$p(x,t)$')
    title(['Probability of particle location at $t=',num2str(tvec(i), '%0.4f'),'$'])
    drawnow;
    % pause(0.005);
end

%% Complex animiation of solution
% close all;
figure
ymax = 1.5*max(abs(soln).^2, [], 'all');
for i=1:tsteps+1
    hold off
    plot3(x, imag(soln(i,:)),real(soln(i,:)), LineWidth=3);
%     scatter3(x, real(soln(i,:)), imag(soln(i,:)),'filled');
    hold on
    % plot3(x(2:N), zeros(length(V)), 0.01*V,'-r',lineWidth=2);
    hold off
    ylim([-1.1*ymax, 1.1*ymax]);
    zlim([-1.1*ymax, 1.1*ymax]);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel('Imaginary', 'Interpreter', 'latex', 'FontSize', 16)
    zlabel('Real', 'Interpreter', 'latex', 'FontSize', 16)
    title(['Wave function at $t=',num2str(tvec(i), '%0.4f'),'$'],'Interpreter','latex', 'FontSize', 16)
    drawnow;
    % exportgraphics(gcf,'new3dcurvewpotential.gif','Append',true);

%     pause(0.005);
end


%%

ymax = max(abs(soln).^2, [], 'all');
figure;
for i = 1:tsteps + 1


    % Plot the probability distribution
    plot(x, abs(soln(i, :)).^2, 'LineWidth', 2);
    hold on
    plot(x(2:N), 0.01 * V, '-r')
    hold off
    ylim([-0.1 * ymax, 1.1 * ymax]);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel('$p(x,t)$', 'Interpreter', 'latex', 'FontSize', 16)
    title(['Probability of particle location at $t=', num2str(tvec(i), '%0.4f'), '$'], 'Interpreter', 'latex', 'FontSize', 16)
    
    exportgraphics(gcf,'testAnimated.gif','Append',true);

end

%%

function V = pot(x)
    V = zeros(1, length(x));
    for i = 1:length(x)
        if (x(i) >= -2) && (x(i) <= 2)
            V(i) = 4;
        elseif x(i) > 2
            V(i) = 1;
        end
    end
end
