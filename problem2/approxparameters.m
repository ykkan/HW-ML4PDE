function [Mf,Mg,Q,c,w,n] = approxparameters(rhomax)
    global T;
    n=1:1:rhomax;
    Q=zeros(rhomax);
    Mf=Q;
    Mg=zeros(rhomax,rhomax+1);
    for rho=1:rhomax
        for k=1:n(rho)
            Q(rho,k)=round(inverse_gamma(rho^(k/2)));
            Mf(rho,k)=round((rho)^((k)/2));
            Mg(rho,k)=rho^(k-1);
        end
        Mg(rho,rho+1)=rho^rho;
    end
    qmax=max(max(Q));
    c=zeros(qmax);
    w=c;
    for k=1:qmax
        [ctemp,wtemp] = lgwt(k,0,T);
        c(:,k)=[flip(ctemp);zeros(qmax-k,1)];
        w(:,k)=[flip(wtemp);zeros(qmax-k,1)];
    end
end