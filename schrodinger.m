%% Solve system

clear; close all;

% Initialize variables
N = 256; L = 10; [D,x] = cheb(N);
x = x*L; D2 = D^2; D2 = D2(2:N,2:N);
[~,w] = clencurt(N); w = w*L;

% Specify initial condition and potential function
psi0 = exp(-(x+5).^2/2).*exp(-1i*1.*x); % IC
% psi0 = exp(-(x+5).^2/2);
% V = 30*exp(-0.5*(x-10).^2) + 30*exp(-0.5*(x+10).^2); % V
% V = 0.5*x.^2;
% V = zeros(N+1);
V = 10*pot(x);

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
dt = 0.001; tf = 0.25; tvec = 0:dt:tf;
tsteps = tf/dt;
soln = zeros(N+1,tsteps+1);
soln(:,1) = psi0;
for i=1:tsteps
    psi = sum(P*(c.*exp(1i*E*tvec(i))), 2);
    soln(:,i+1) = [0;psi;0];
end
soln = soln';

%% Plot first couple eigenfunctions (numerical)

figure
hold on
for i=1:4
    plot(x,[0;P(:,i);0],LineWidth=2);
end
legend('$n=1$','$n=2$','$n=3$','$n=4$',Interpreter='latex')

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

%% Animiation of solution

figure
ymax = max(abs(soln).^2, [], 'all');
for i=1:tsteps+1
    plot(x, abs(soln(i,:)).^2, LineWidth=2);
    hold on
    plot(x(2:N), 0.01*V, '-r')
    hold off
    ylim([-0.1*ymax, 1.1*ymax]);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel('$p(x,t)$', 'Interpreter', 'latex', 'FontSize', 16)
    title(['Probability of particle location at $t=',num2str(tvec(i), '%0.4f'),'$'],'Interpreter','latex', 'FontSize', 16)
    drawnow;
%     pause(0.005);
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
    plot3(x(2:N), zeros(length(V)), 0.01*V,'-r',lineWidth=2);
    hold off
    ylim([-1.1*ymax, 1.1*ymax]);
    zlim([-1.1*ymax, 1.1*ymax]);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16)
    ylabel('Imaginary', 'Interpreter', 'latex', 'FontSize', 16)
    zlabel('Real', 'Interpreter', 'latex', 'FontSize', 16)
    title(['Wave function at $t=',num2str(tvec(i), '%0.4f'),'$'],'Interpreter','latex', 'FontSize', 16)
%     drawnow;
    exportgraphics(gcf,'new3dcurvewpotential.gif','Append',true);

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
