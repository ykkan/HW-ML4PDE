clear all;
u = @(x,y) sin(pi*x)*sin(pi*y);
f = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);
L1 = @(A,B,h) h*sum(abs(A-B));
L2 = @(A,B,h) sqrt(sum((A-B).^2*h));
Linf = @(A,B) max(abs(A-B));


ns = [51, 101, 201, 401, 601, 801];
errors_L1 = zeros(length(ns),1);
errors_L2 = zeros(length(ns),1);
errors_Linf = zeros(length(ns),1);
hs = zeros(length(ns),1);
for iter=1:length(ns)
   n=ns(iter);
   U=poisson2d(f,n);
   U_exact=zeros(n^2,1);
   x=linspace(0,1,n);
   for j=1:n
       for i=1:n
           k=i+(j-1)*n;
           U_exact(k)=u(x(i),x(j));
       end
   end
   h=1/(n-1);
   errors_L1(iter)=L1(U,U_exact,h);
   errors_L2(iter)=L2(U,U_exact,h);
   errors_Linf(iter)=Linf(U,U_exact);
   hs(iter)=h;
end

xx=log(hs);
yy=log(errors_L2);
P=polyfit(xx,yy,1);
g = @(x) exp(P(2)).*x.^P(1);

loglog(hs, errors_L2,'o',hs,g(hs), '-')
xlabel("h",'FontSize',20)
ylabel("l^2 error",'FontSize',20)
legend({'data','error(h) = 0.8738\timesh^{0.5010}'},'Location','northwest','FontSize',20)
set(gcf, 'PaperPosition', [0 0 10.66 6]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
set(gcf, 'PaperSize', [10.66 6]); %Keep the same paper size
saveas(gcf, 'convergence', 'pdf')