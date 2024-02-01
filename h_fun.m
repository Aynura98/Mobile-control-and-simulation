function[y]=h_fun(x,u)      % Get the current distances robot/beacons
parameters;
P=x(1:2)';      % Get current coordinates of the robot
y=zeros(size(S,1),1);   % Preallocation of variable
    for i=1:length(y)
        y(i)=norm(S(i,:)-P);    % Find the distance between the robot position and each one of the beacons
    end
end