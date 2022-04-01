function y=inverse_gamma(x)
    c=0.036534;
    L= log((x+c)/sqrt(2*pi));
    y=L/lambertw(L/exp(1))+0.5;
end