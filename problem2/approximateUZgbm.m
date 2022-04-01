function [u, z] = approximateUZgbm(n,rho,x,s)
 global Mf Mg Q c w T dim f g mu sigma;
 cloc=(T-s)*c/T+s;
 wloc=(T-s)*w/T;
 MC=Mg(rho,n+1);
 W=sqrt(T-s)*randn(dim,MC);
 X=repmat(x,1,MC).*exp((mu-sigma^2/2)*(T-s)+sigma*W);
 xi=g(X);
 u=sum(xi,2)/MC;
 z=sum(repmat(xi-g(x)*ones(1,MC),dim,1).*W,2)./(MC*(T-s));
 for l=0:(n-1)
  q=Q(rho,n-l);
  d=cloc(1:q,q)-[s;cloc(1:(q-1),q)];
  MC=Mf(rho,n-l);
  X=repmat(x,1,MC);
  W=zeros(dim,MC);
  for k=1:q
   dW=sqrt(d(k))*randn(dim,MC);
   W=W+dW;
   X=X.*exp((mu-sigma^2/2)*d(k)+sigma*dW);
   [U, Z]=cellfun(@approximateUZgbm, num2cell(l*ones(1,MC)), num2cell(rho*ones(1,MC)),...
       num2cell(X,1), num2cell(cloc(k,q)*ones(1,MC)),'UniformOutput',false);
   y=f(cloc(k,q),X,cell2mat(U),cell2mat(Z));
   u=u+wloc(k,q)*sum(y)/MC;
   z=z+wloc(k,q).*sum(repmat(y,dim,1).*W,2)./(MC*(cloc(k,q)-s));
    if l>0
     [U, Z]=cellfun(@approximateUZgbm, num2cell((l-1)*ones(1,MC)), num2cell(rho*ones(1,MC)),...
         num2cell(X,1), num2cell(cloc(k,q)*ones(1,MC)),'UniformOutput',false);
     y=f(cloc(k,q),X,cell2mat(U),cell2mat(Z));
     u=u-wloc(k,q)*sum(y)/MC;
     z=z-wloc(k,q).*sum(repmat(y,dim,1).*W,2)./(MC*(cloc(k,q)-s));
    end
  end  
 end
end