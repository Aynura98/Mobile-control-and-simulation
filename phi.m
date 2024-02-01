function [Xk1]=phi(Xk,uk)
parameters;

xk=Xk(1);
yk=Xk(2);
thetak=Xk(3);
vk=uk(1);
wk=uk(2);

Xk1=[xk+(vk*Ts*cos(thetak+wk*Ts));
    yk+(vk*Ts*sin(thetak+wk*Ts));
    thetak+wk*Ts
    ];
end