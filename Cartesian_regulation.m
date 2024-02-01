clear;
close all;
clc;
%% Reference trajectory
Path_planning;      

load Path_planning_variables.mat X dX nsteps O T;    

parameters;

X = X';  
dX = dX';
X0 = [X(1,1),X(2,1),X(3,1)];  
x_star = X(1,1);   
y_star = X(2,1);
theta_star = X(3,1); 
xd_star = dX(1,1);   
yd_star = dX(2,1);

k1 = 5;
k2 = 10;

X0_column_size = size(X0,2);
step_size = nsteps-1;
XReal = zeros(step_size,X0_column_size);          % state  
X_Kestimated = zeros(X0_column_size,step_size);   % Kalman estimator   
X_Kpredredicted = zeros(X0_column_size,step_size);% Kalman predictor   
v = zeros(step_size,1);                           % translational speed input
w = zeros(step_size,1);                           % rotational speed input

v_noise = zeros(nsteps,size(V,1));    
w_noise = zeros(nsteps,size(W,1));    
for k = 1:nsteps       
    w_noise(k,:) = (sqrt(W) * randn(size(W,1),1))';   % proccess noise
    v_noise(k,:) = (sqrt(V) * randn(size(V,1),1))';   % measurement noise
end

%% Robot evolution
h = figure(6);          

Xhat0 = X0';

ekf = EKF(@phi,@h_fun,@Ak_LZ,@Ck_LZ,Xhat0,P0,V,W); % object

f_robot=@(t,X)cartesian_reg(X,x_star,y_star,theta_star,k1,k2);
[t,evolution] = ode45(f_robot,[0,Ts],X0);

for i = 1:length(X)-1
    evolution_last_row = evolution(length(evolution),:);
    XReal(i,:) = evolution_last_row + w_noise(i,:); 
    %XReal(i,:) = evolution(length(evolution),:); 
    thetareal = enrollTheta(XReal(i,3));    % to keep theta between 0 and 2pi
    state_real = [XReal(i,1:2) thetareal];
    Xdot=cartesian_reg([XReal(i,1:2) thetareal],x_star,y_star,theta_star,k1,k2);
    v(i,1) = norm(Xdot(1,1) + Xdot(2,1));   %translational speed input          
    w(i,1) = Xdot(3,1);                     %rotational speed input
    distance = h_fun([XReal(i,1);XReal(i,2)],0);
    uk = [v(i,1) w(i,1)];                   %controller output
    yk = distance + v_noise(i,:)';          %measurement
    ekf.update(yk,uk);                      %update
    [Xhat_k1k,Pk1k] = ekf.get_k1k_data();   %obtain prediction result
    [Xhat_kk,Pkk] = ekf.get_kk_data();      %obtain estimation result 
    X_Kestimated(:,i) = Xhat_kk(:,i);  
    X_Kpredredicted(:,i) = Xhat_k1k(:,i);
    Xhatk1k = X_Kpredredicted(:,i);
    x_star = X(1,i + 1);                    %next points   
    y_star = X(2,i + 1);        
    xd_star = dX(1,i + 1);
    yd_star = dX(2,i + 1);
    
    f_robot=@(t,X)cartesian_reg(X,x_star,y_star,theta_star,k1,k2);
    [t,evolution] = ode45(f_robot,[0,Ts],X_Kestimated(:,i));
    
    % Plotting simulation
    figure(h);
    cla(h);     
    hold on
    grid on;
    scatter(S(:,1),S(:,2),50,"square",'MarkerFaceColor','k','MarkerEdgeColor','k'); % beacons
    plot(X(1,:),X(2,:),'--k','LineWidth',1);    % reference trajectory - dashed black line 
    
    plot(O(:,1),O(:,2),'ok','LineWidth',2);     % obstacles
    plot(X(size(X,1),1),X(size(X,1),2),'xb','LineWidth',2); % goal point
    
    plot_sensors_data(XReal(i,:),yk,h);     % beacons - dashed red lines  
    plot_robot_traj(XReal',i,1,[],h);       % real evolution - blue line
    plot_robot_traj(X_Kestimated,i,2,[],h); % estimated evolution - green line
    
end

%% Results
figure('Name','Trajectories');
for i = 1:3
    subplot(3,1,i);
    hold on;
    grid on
    title(['State component ', num2str(i)]);
    plot(X(i,:),'--k','LineWidth',1);
    plot(XReal(:,i),'r','LineWidth',2);
    plot(X_Kpredredicted(i,:),'--g','LineWidth',1);
    plot(X_Kestimated(i,:),'--b','LineWidth',2);
    legend('Reference','Real','Predicted','Estimated','Location','best');
end

figure('Name','Error of estimation and prediction');
subplot(3,1,1);
title('Estimation error')
hold on;
grid on
plot(X_Kestimated(1,:)-X(1,1:end-1),'--b','LineWidth',1);
plot(X_Kestimated(2,:)-X(2,1:end-1),'--r','LineWidth',1);
plot(zeros(size(X_Kpredredicted,2)),'-k');  % To have a reference line
legend('First component','Second component','Third component','Location','best');

subplot(3,1,2);
title('Prediction error')
hold on;
grid on
plot(X_Kpredredicted(1,:)-X(1,1:end-1),'--b','LineWidth',1);
plot(X_Kpredredicted(2,:)-X(2,1:end-1),'--r','LineWidth',1);
plot(zeros(size(X_Kpredredicted,2)),'-k');
legend('First component','Second component','Third component','Location','best');

subplot(3,1,3);
title('Estimation vs Prediction')
hold on;
grid on
plot(X_Kpredredicted(1,:)-X_Kestimated(1,:),'--b','LineWidth',1);
plot(X_Kpredredicted(2,:)-X_Kestimated(2,:),'--r','LineWidth',1);
plot(zeros(size(X_Kpredredicted,2)),'-k');
legend('First component','Second component','Third component','Location','best')
hold off

%% Enroll Theta
function[theta] = enrollTheta(theta)
    for i=1:length(theta)
        theta(i) = mod(theta(i),2 * pi); 
        if theta(i)>pi
            theta(i) = theta(i)-2 * pi;
        end
    end
end

%% Plot robot trajectory
function[h] = plot_robot_traj(Xtraj,k,robot_type,axis_vec,h)
    if nargin < 5 || isempty(h)
        h = figure();
    end
    if nargin < 4
        axis_vec = [];
    end
    
    if robot_type == 1
        cc = 'b';
    else
        cc = '--g';
    end
    if k > 0
        x = Xtraj(:,1:k);
        hold on;
        plot(x(1,:),x(2,:),cc,'LineWidth',2);
        plot_robot(x(:,end),robot_type,h);
    else
        for kk = 1:size(Xtraj,2)
            cla;
            plot_robot_traj(Xtraj,kk,robot_type,axis_vec,h);
            pause(0.05);
        end
    end
    if isempty(axis_vec) == 0
        axis(axis_vec);
    end
end

%% Plot robot
function[h] = plot_robot(Xk,robot_type,h_fig)
    d = 0.5;
    if nargin < 3
        h = figure();
    else
        h = h_fig;
        figure(h);
    end
    
    if robot_type == 1
        color1 = [0 0 1];
        color2 = [1 0 0];
    elseif robot_type == 2
        color1 = [0 1 0];
        color2 = [1 0 0];
    else
        color1 = [1 0 1];
        color2 = [1 0 0];
    end
    
    P1 = [Xk(1) + d/2 * cos(pi/2 + Xk(3)) Xk(2) + d/2 * sin(pi/2 + Xk(3))];
    P2 = [Xk(1) + d/2 * cos(-pi/2 + Xk(3)) Xk(2) + d/2 * sin(-pi/2 + Xk(3))];
    P3 = [Xk(1) + d * cos(Xk(3)) Xk(2) + d * sin(Xk(3))];
    
    plot(Xk(1),Xk(2),'o','color',color2,'LineWidth',2);
    hold on;
    line([P1(1) P3(1)],[P1(2) P3(2)],'color',color1,'LineWidth',2);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color',color1,'LineWidth',2);
    line([P2(1) P3(1)],[P2(2) P3(2)],'color',color1,'LineWidth',2);
end

%%
function[h] = plot_sensors_data(X,y,h)
    if nargin == 2 || isempty(h)
        h = figure();
    end
    parameters;
    hold on;
    for i = 1:length(y)
        theta = atan2(X(2) - S(i,2),X(1) - S(i,1));
        plot([S(i,1) S(i,1) + y(i) * cos(theta)],[S(i,2) S(i,2) + y(i) * sin(theta)],'--r','LineWidth',2);
    end
end
%%
function Xdot = cartesian_reg(X, x_star,y_star,theta_star,k1,k2)
x = X(1);
y = X(2);
theta = X(3);

ex = x_star-x;
ey = y_star-y;

v = k1*(ex*cos(theta)+ey*sin(theta));
% w=k2*(atan2(ey,ex)-theta);
if norm([ex;ey]) <= 0.001
    w = 0;
else
    w = k2 * (delta_angle(atan2(ey,ex),theta));
end


Xdot = [v * cos(theta);
      v * sin(theta);
      w];

end

