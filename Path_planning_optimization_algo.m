%%
clear;
close all;
clc;
%%
vortex=0;
% Initial position
X0 = [1;1];
% Goal position
G = [14;14];
ob1_x = 0:0.5:2;
ob1_y = 2*ones(1,size(ob1_x,2));
ob2_x = 5:0.5:8;
ob2_y = 7*ones(1,size(ob2_x,2));
ob3_x = 2:0.5:4;
ob3_y = 0*ones(1,size(ob3_x,2));
ob4_y = 2:0.5:4;
ob4_x = 8*ones(1,size(ob4_y,2));
ob5_x = 9:0.5:12;
ob5_y = 11*ones(1,size(ob5_x,2));
ob6_y = 9:0.5:12;
ob6_x = 15*ones(1,size(ob6_y,2));
O = vertcat([ob1_x' ob1_y'],[ob2_x' ob2_y'],[ob3_x' ob3_y'],[ob4_x' ob4_y'],[ob5_x' ob5_y'],[ob6_x' ob6_y']) ;


rotations=ones(size(O,1),1);

dmin=1.2;

Ja=@(x,y,Gx,Gy)((1/2)*((x-Gx).^2+(y-Gy).^2));% potential Ja
Jr=@(x,y,Ox,Oy)(1./((x-Ox).^2+(y-Oy).^2));% potential Jr


nablaJaX=@(x,y,Gx,Gy)(x-Gx);% gradient of Ja w.r.t. X
nablaJaY=@(x,y,Gx,Gy)(y-Gy);% gradient of Ja w.r.t. Y
nablaJrX=@(x,y,Ox,Oy)(2*(Ox-x)./(((x-Ox).^2+(y-Oy).^2).^2));% derivative of Jr w.r.t. X
nablaJrY=@(x,y,Ox,Oy)(2*(Oy-y)./(((x-Ox).^2+(y-Oy).^2).^2));% derivative of Jr w.r.t. Y


rho=@(x,y,Ox,Oy)((x-Ox).^2+(y-Oy).^2<=dmin^2);% defines field where repulsion ...
% effect can be detected, in this way we can avoid the minimum position of overall...
% potential is changed due to repulsive potential
% returns 1 if the distance between defined obstacle and position is less
% than the given threshold
% returns 0, otherwise

wa=2;
wo=5;
Tend=130;
deltaXY=0.2;
%%
xm=min(X0(1),G(1));xm=min(xm,min(O(:,1)));
xM=max(X0(1),G(1));xM=max(xM,max(O(:,1)));
ym=min(X0(2),G(2));ym=min(ym,min(O(:,2)));
yM=max(X0(2),G(2));yM=max(yM,max(O(:,2)));
%%
figure();hold on;grid on;
plot(X0(1),X0(2),'xm','LineWidth',2);
plot(G(1),G(2),'xr','LineWidth',2);
plot(O(:,1),O(:,2),'ok','LineWidth',2);
axis ([0 14 0 14]);
axis('equal');
%%
xx=xm-2:deltaXY:xM+2;
yy=ym-2:deltaXY:yM+2;
[XX,YY] = meshgrid(xx,yy);

Za=Ja(XX,YY,G(1),G(2));
nablaJaXX=nablaJaX(XX,YY,G(1),G(2));
nablaJaYY=nablaJaY(XX,YY,G(1),G(2));

figure('Name','Attraction');
subplot(121);hold on;grid on;surf(XX,YY,Za);
subplot(122);hold on;grid on;quiver(XX,YY,-nablaJaXX,-nablaJaYY);

Zr=zeros(size(Za));
nablaJrXX=zeros(size(nablaJaXX));
nablaJrYY=zeros(size(nablaJaXX));

for i=1:size(O,1)
    oi=O(i,:);
    Zr=Zr+Jr(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2));% to ensure piecewise definition

    if ~vortex
        nablaJrXX=nablaJrXX+nablaJrX(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2))*rotations(i);
        nablaJrYY=nablaJrYY+nablaJrY(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2))*rotations(i);
    else
        nablaJrXX=nablaJrXX-nablaJrY(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2))*rotations(i);
        nablaJrYY=nablaJrYY+nablaJrX(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2))*rotations(i);
    end
end

figure('Name','Repulsion');
subplot(121);hold on;grid on;surf(XX,YY,Zr);
subplot(122);hold on;grid on;quiver(XX,YY,-nablaJrXX,-nablaJrYY);

J=wa*Za+wo*Zr;
nablaJx=wa*nablaJaXX+wo*nablaJrXX;
nablaJy=wa*nablaJaYY+wo*nablaJrYY;

nablaJxn=nablaJx./sqrt(nablaJx.^2+nablaJy.^2);% normalized version in order to show the directions of arrows
nablaJyn=nablaJy./sqrt(nablaJx.^2+nablaJy.^2);


figure('Name','Overall');
subplot(2,2,[1 3]);hold on;grid on;surf(XX,YY,J);
subplot(222);hold on;grid on;quiver(XX,YY,-nablaJx,-nablaJy);
subplot(224);hold on;grid on;quiver(XX,YY,-nablaJxn,-nablaJyn);
%%

figure();
hold on;grid on;
quiver(XX,YY,-nablaJxn,-nablaJyn);
plot(X0(1),X0(2),'xm','LineWidth',2);
plot(G(1),G(2),'xr','LineWidth',2);
plot(O(:,1),O(:,2),'ok','LineWidth',2);

%% APF
alpha=0.02;
% alpha=0.1;
TH=0.2;
nIter=100;
X=zeros(2,nIter);
X(:,1)=X0;
for k=2:nIter
    Xcur=X(:,k-1);
    nablaA=[
            nablaJaX(Xcur(1),Xcur(2),G(1),G(2));
            nablaJaY(Xcur(1),Xcur(2),G(1),G(2))
        ];
    nablaR=[0;0];
    for i=1:size(O,1)
        oi=O(i,:);
        if ~vortex
            nablaR=nablaR+[
                nablaJrX(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),Xcur(2),oi(1),oi(2))*rotations(i);
                nablaJrY(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),Xcur(2),oi(1),oi(2))*rotations(i)
                ];
        else
            nablaR=nablaR+[0 -1;1 0]*[
                nablaJrX(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),Xcur(2),oi(1),oi(2))*rotations(i);
                nablaJrY(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),Xcur(2),oi(1),oi(2))*rotations(i)
                ];
        end
    end
    nablaJ=wa*nablaA+wo*nablaR;
    Xsucc=Xcur-alpha*nablaJ;
    X(:,k)=Xsucc;
    if norm(Xsucc-G)<=TH
        break;
    end
end
if k<nIter
    X(:,k+1:end)=[];
end

for i=1:size(X,2)
    plot(X(1,i),X(2,i),'xk','LineWidth',2);
end
