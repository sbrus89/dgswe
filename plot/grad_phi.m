function gp = grad_phi(px,h)

gp = (-8*px*(1-sum(px.*px)/h^2)^3)/h^2;
