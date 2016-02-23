function xg = gaussian_function(x,sigma)

xg = 1/(sigma*sqrt(2*pi))*exp(-x.^2/(2*(sigma^2)));

