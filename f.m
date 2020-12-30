function func = f(x)

global delta;
% discrete f(x), LHS of DE
funcdef = @(x) [(abs(x)<=0.1) (abs(0.2-abs(x))<=0.1) (abs(0.4-abs(x))<=0.1) (abs(0.6-abs(x))<=0.1) (abs(0.9-abs(x))<=0.2) (abs(x)>=1.1) ]*(1/delta)*[-5; 25; 1; -30; 20; 1];

func = funcdef(x);
end
