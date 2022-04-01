function y=approximateUfinitediffabm(t0,x0,N)
    [T,dim,f,g,eta,mu,sigma]=modelparameters();
    h=(T-t0)/N;
    t=t0:h:T;
    u=mu*h+sigma*sqrt(h);
    d=mu*h-sigma*sqrt(h);
    x=x0+d*N+(0:N)*(u-d);
    M=(1/2*[full(gallery('tridiag',ones(N-1,1),ones(N,1),zeros(N-1,1)));[zeros(1,N-1),1]]);
    L=1/(2*sqrt(h))*([full(gallery('tridiag',ones(N-1,1),-ones(N,1),zeros(N-1,1)));[zeros(1,N-1),1]]);
    y=g(x);
    for i=N:-1:1
        x=x(1:i)-d;
        z=y*L(1:i+1,1:i);
        y=y*M(1:i+1,1:i);
        y=y+h*f(t(i),x,y,z);
    end
end