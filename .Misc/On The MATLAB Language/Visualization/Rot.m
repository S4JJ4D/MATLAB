function R = Rot(axis,angle)
%% This function receives rotation axis and angle of rotation(expressed in degrees) as input and returns their corresponding "Rotation Matrix".
% x = deg2rad(angle); %angle MUST be in degrees.
% axis = upper(axis);
x = angle;

if strcmpi(axis,'X')
    
    R = [1 0 0; 0 cos(x) -sin(x) ; 0 sin(x) cos(x)];
    
elseif strcmpi(axis,'Y')
    
    R = [cos(x) 0 sin(x) ; 0 1 0 ; -sin(x) 0 cos(x)];
    
elseif strcmpi(axis,'Z')
    
    R = [cos(x) -sin(x) 0 ; sin(x) cos(x) 0 ; 0 0 1];
    
else 
    
    error('Invalid Input');
    
end


end