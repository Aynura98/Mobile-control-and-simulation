function[Ak]=Ak_LZ(X_hat,u)
parameters;

vk=u(1);
wk=u(2);
theta=X_hat(3);

Ak=[
    1 0 -vk*Ts*sin(theta+wk*Ts);
    0 1  vk*Ts*cos(theta+wk*Ts);
    0 0             1
    ];
end