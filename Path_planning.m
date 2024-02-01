clear;
close all;
clc;

vortex = 0; 
T = 130;   
% Initial position
X0 = [1;1];
% Goal position
G = [14;14];

%% Obstacles position
% Addition of walls
ob1_x = 0:0.5:2;
ob1_y = 2 * ones(1,size(ob1_x,2));
ob2_x = 5:0.5:8;
ob2_y = 7 * ones(1,size(ob2_x,2));
ob3_x = 2:0.5:4;
ob3_y = 0.5 * ones(1,size(ob3_x,2));
ob4_y = 2:0.5:5;
ob4_x = 7 * ones(1,size(ob4_y,2));
ob5_x = 9:0.5:12;
ob5_y = 11 * ones(1,size(ob5_x,2));

O = vertcat([ob1_x' ob1_y'],[ob2_x' ob2_y'],[ob3_x' ob3_y'],[ob4_x' ob4_y'],[ob5_x' ob5_y']) ;

% Direction of rotation
rotations = ones(size(O,1),1);      

dmin = 1;

Ja = @(x,y,Gx,Gy)((1/2) * ((x - Gx).^2 + (y - Gy).^2));% potential Ja
Jr = @(x,y,Ox,Oy)(1./((x - Ox).^2 + (y-Oy).^2));% potential Jr

nablaJaX = @(x,y,Gx,Gy)(x - Gx);% derivative of Ja w.r.t. X
nablaJaY = @(x,y,Gx,Gy)(y - Gy);% derivative of Ja w.r.t. Y
nablaJrX = @(x,y,Ox,Oy)(2 * (Ox - x)./(((x - Ox).^2+(y - Oy).^2).^2));% derivative of Jr w.r.t. X
nablaJrY = @(x,y,Ox,Oy)(2 * (Oy - y)./(((x - Ox).^2+(y - Oy).^2).^2));% derivative of Jr w.r.t. Y

rho = @(x,y,Ox,Oy)((x - Ox).^2 + (y - Oy).^2 <= dmin^2);

alpha = 0.01;
deltaXY = 0.2;    
                      
wa = 5; % reaches goal faster passing obstacles in close distance                           
wo = 2;  
%%
xm = min(X0(1),G(1)); xm=min(xm,min(O(:,1)));
xM = max(X0(1),G(1)); xM=max(xM,max(O(:,1)));
ym = min(X0(2),G(2)); ym=min(ym,min(O(:,2)));
yM = max(X0(2),G(2)); yM=max(yM,max(O(:,2)));

figure('Name','Trajectory');
hold on;
grid on;
plot(X0(1),X0(2),'xr','LineWidth',2);       % Starting point
plot(G(1),G(2),'xb','LineWidth',2);         % Goal point
plot(O(:,1),O(:,2),'ok','LineWidth',2);     % Obstacles

%%
xx = xm - 2:deltaXY:xM + 2;
yy = ym - 2:deltaXY:yM + 2;
[XX,YY] = meshgrid(xx,yy);

Za = Ja(XX,YY,G(1),G(2));
nablaJaXX = nablaJaX(XX,YY,G(1),G(2));
nablaJaYY = nablaJaY(XX,YY,G(1),G(2));

figure('Name','Attraction');
subplot(121); hold on; grid on; surf(XX,YY,Za);
subplot(122); hold on; grid on; quiver(XX,YY,-nablaJaXX,-nablaJaYY);

Zr = zeros(size(Za));
nablaJrXX = zeros(size(nablaJaXX));
nablaJrYY = zeros(size(nablaJaXX));

for i = 1:size(O,1)
    oi = O(i,:);
    Zr = Zr + Jr(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2));% to ensure piecewise definition

    if ~vortex
        nablaJrXX = nablaJrXX + nablaJrX(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2)) * rotations(i);
        nablaJrYY = nablaJrYY + nablaJrY(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2)) * rotations(i);
    else
        nablaJrXX = nablaJrXX - nablaJrY(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2)) * rotations(i);
        nablaJrYY = nablaJrYY + nablaJrX(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2)) * rotations(i);
    end
end

figure('Name','Repulsion');
subplot(121);hold on;grid on;surf(XX,YY,Zr);
subplot(122);hold on;grid on;quiver(XX,YY,-nablaJrXX,-nablaJrYY);

J = wa * Za + wo * Zr;
nablaJx = wa * nablaJaXX + wo * nablaJrXX;
nablaJy = wa * nablaJaYY + wo * nablaJrYY;

nablaJxn = nablaJx./sqrt(nablaJx.^2 + nablaJy.^2);% normalized version in order to show the directions of arrows
nablaJyn = nablaJy./sqrt(nablaJx.^2 + nablaJy.^2);

figure('Name','Overall');
subplot(2,2,[1 3]);hold on;grid on;surf(XX,YY,J);
subplot(222);hold on;grid on;quiver(XX,YY,-nablaJx,-nablaJy);
subplot(224);hold on;grid on;quiver(XX,YY,-nablaJxn,-nablaJyn);

%%
% Calculate trajectory
[t,X] = ode45(@(t,X)(nablaEvolution(t,X,nablaJaX,nablaJaY,nablaJrX,nablaJrY,rotations,vortex,wa,wo,G,O,alpha,rho)),[0,T],X0);
nsteps = size(X,1);
dX = zeros(nsteps,2);

x_coordinates = X(:,1);
y_coordinates = X(:,2);
for i = 1:size(t,1)     
    Xdot = nablaEvolution(t(i),[x_coordinates(i),y_coordinates(i)],nablaJaX,nablaJaY,nablaJrX,nablaJrY,rotations,vortex,wa,wo,G,O,alpha,rho);
    dX(i,:) = [Xdot(1) Xdot(2)];
end

if (nsteps <= 2 * T)           
    X_ref = zeros(nsteps * 2,size(X,2));
    dX_ref = zeros(nsteps * 2,size(dX,2));
    for i = 1:size(dX,2)        
        X_ref(:,i)   = interp(X(:,i),2);
        dX_ref(:,i)  = interp(dX(:,i),2);
    end
    nsteps = size(X_ref,1);    
    theta  = atan2(dX_ref(:,2),dX_ref(:,1));
    X = [X_ref theta];
    dX = dX_ref;
elseif (nsteps >= 3.75 * T)    
    X_ref = zeros(ceil(nsteps/2),size(X,2));
    dX_ref = zeros(ceil(nsteps/2),size(dX,2));
    for i = 1:size(dX,2)        
        X_ref(:,i) = decimate(X(:,i),2);
        dX_ref(:,i) = decimate(dX(:,i),2);
    end
    nsteps = size(X_ref,1);    
    theta = atan2(dX_ref(:,2),dX_ref(:,1));
    X = [X_ref theta];
    dX = dX_ref;
end

figure("Name","Trajectory overall potential fields");
hold on;grid on;
quiver(XX,YY,-nablaJxn,-nablaJyn);
plot(X0(1),X0(2),'xm','LineWidth',2);
plot(G(1),G(2),'xr','LineWidth',2);
plot(O(:,1),O(:,2),'ok','LineWidth',2);
plot(X(:,1),X(:,2),'k','LineWidth',2);

save Path_planning_variables X dX nsteps O T 
clearvars   j

% Function implementation
function[dX] = nablaEvolution(~,X,nablaJaX,nablaJaY,nablaJrX,nablaJrY,rotations,vortex,wa,wo,G,O,alpha,rho)
    x = X(1);   y = X(2);
    nablaJx = wa * nablaJaX(x,y,G(1),G(2));
    nablaJy = wa * nablaJaY(x,y,G(1),G(2));
    for k = 1:size(O,1)
        if ~vortex
            nablaJx = nablaJx + wo * nablaJrX(x,y,O(k,1),O(k,2)) * rho(x,y,O(k,1),O(k,2)) * rotations(k);
            nablaJy = nablaJy + wo * nablaJrY(x,y,O(k,1),O(k,2)) * rho(x,y,O(k,1),O(k,2)) * rotations(k);            
        else 
            nablaJr = [nablaJrX(x,y,O(k,1),O(k,2)) * rotations(k);
                     nablaJrY(x,y,O(k,1),O(k,2)) * rotations(k)];
            nalbaJr = [1 0;0 -1] * nablaJr;
            nablaJx = nablaJx + wo * nalbaJr(1);
            nablaJy = nablaJy + wo * nalbaJr(2);
        end

    end
    dX = -alpha * [nablaJx;nablaJy];
end