function [t0, x0, u0] = shift(T, t0, x0, u,f)
% add noise to the control actions before applying it
con_cov = diag([0.005 deg2rad(2)]).^2;
con = u(1,:)' + sqrt(con_cov)*randn(2,1); 
st = x0;

f_value = f(st,con);   
st = st+ (T*f_value);

x0 = full(st);
t0 = t0 + T;

u0 = [u(2:size(u,1),:);u(size(u,1),:)]; % shift the control action 
end