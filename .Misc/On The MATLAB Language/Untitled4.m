%% 
hws = get_param('Robotics_HW8','modelworkspace');
value = str2double(get_param(gcbh,'Sub_theta_1'));

hws.assignin('theta_1',value);



%% 

p = Simulink.Mask.get(gcbh);
hws = get_param('Robotics_HW8','modelworkspace');
value = get_param(gcbh,'InvActivationCheck');

if strcmp(value,'on')
    hws.assignin('isActivated',1);
    p.getParameter('X_ii').set('Enabled','on')
    p.getParameter('Y_ii').set('Enabled','on')
    p.getParameter('Z_ii').set('Enabled','on')
    p.getParameter('X_ff').set('Enabled','on')
    p.getParameter('Y_ff').set('Enabled','on')
    p.getParameter('Z_ff').set('Enabled','on')
else
    hws.assignin('isActivated',0);
    p.getParameter('X_ii').set('Enabled','off')
    p.getParameter('Y_ii').set('Enabled','off')
    p.getParameter('Z_ii').set('Enabled','off')
    p.getParameter('X_ff').set('Enabled','off')
    p.getParameter('Y_ff').set('Enabled','off')
    p.getParameter('Z_ff').set('Enabled','off')
end
%hws.assignin('X_i',value);

%% hws = get_param('Robotics_HW8','modelworkspace');
value = str2double(get_param(gcbh,'X_ff'));

hws.assignin('X_f',value);
