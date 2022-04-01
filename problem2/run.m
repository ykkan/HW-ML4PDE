clear all;
global Mf Mg Q c w T dim f g mu sigma eta;
rng(2017)
format long
average=10;
rhomin=1;
rhomax=6;
[T,dim,f,g,eta,mu,sigma]=modelparameters();
[Mf,Mg,Q,c,w,n] = approxparameters(rhomax);
value=zeros(average,rhomax-rhomin+1);
time=value;
for rho=rhomin:rhomax
  for k=1:average
    tic
    [a, b]=approximateUZgbm(n(rho),rho,100*ones(dim,1),0);
    value(k,rho-rhomin+1)=a;
    time(k,rho-rhomin+1)=toc;
  end
end
name = [datestr(now, 'yymmddTHHMMSS') '.mat']; 
save(name,'n','Q','Mf','Mg','value','time')
plotruntimevsdim(rhomin, rhomax, time)