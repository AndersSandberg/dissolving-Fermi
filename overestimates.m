[V,rho1]=meshgrid(0:0.01:1,1:0.1:10);
surf(V,rho1, (rho1.^2.*V+1-V)./(rho1.*V+(1-V)).^2)
xlabel('V')
ylabel('\rho_1/\rho_2')
zlabel('$\hat{\rho}/\bar{\rho}-1$','Interpreter','latex')
% print -depsc overestimate