function[Ck]=Ck_LZ(X,u)
parameters;
x=X(1);
y=X(2);
theta=X(3);

Ck=zeros(size(S,1),3);
for i=1:size(Ck,1)
    Sx=S(i,1);
    Sy=S(i,2);
    dx=(-1).*(Sx+(-1).*x).*((Sx+(-1).*x).^2+(Sy+(-1).*y).^2).^(-1/2);
    dy=(-1).*((Sx+(-1).*x).^2+(Sy+(-1).*y).^2).^(-1/2).*(Sy+(-1).*y);
    dtheta=0;
    Ck(i,:)=[dx dy dtheta];
end
end