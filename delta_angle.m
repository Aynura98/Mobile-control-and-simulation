function[dAngle]=delta_angle(theta1,theta2)
dAngle=atan2(sin(theta1-theta2),cos(theta1-theta2));
end