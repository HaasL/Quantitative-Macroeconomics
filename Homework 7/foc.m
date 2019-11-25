function fval = foc(cons,x,survival)

global betta tetta r g gridx vpfun epsi probepsi ne nx 

vpp1 = evalvp(cons,x);
margu = margutil(cons);
fval = margu - survival*(betta * (1+g)^(-tetta) * (1+r) * vpp1);

end
