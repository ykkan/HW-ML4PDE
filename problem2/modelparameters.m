function [T,dim,f,g,eta,mu,sigma] = modelparameters()
    T = 0.5;
    dim=100;
    mu=0;
    sigma=0.25;
    f = @(t,x,y,z) sigma*(y-(2+sigma^2*dim)/(2*sigma^2*dim)).*sum(z,1); 
    g = @(x) 1-1./(1+exp(T+sum(x,1)));
    eta=@(x) x;
end