modelWorkspaceHandle = get_param('Robotics_HW8','modelworkspace');

if getVariable(modelWorkspaceHandle,'isActivated')

xi = getVariable(modelWorkspaceHandle,'X_i');
yi = getVariable(modelWorkspaceHandle,'Y_i');
zi = getVariable(modelWorkspaceHandle,'Z_i');
xf = getVariable(modelWorkspaceHandle,'X_f');
yf = getVariable(modelWorkspaceHandle,'Y_f');
zf = getVariable(modelWorkspaceHandle,'Z_f');

[t1,t2,dd] = Inv_Kin_Scara(xi,yi,zi);
[T1,T2,DD] = Inv_Kin_Scara(xf,yf,zf);

ICBlock_handle = Simulink.Mask.get('Robotics_HW8/IC Block');
t1h = ICBlock_handle.getParameter('Sub_theta_1');
t2h = ICBlock_handle.getParameter('Sub_theta_2');
dh = ICBlock_handle.getParameter('Sub_d');

set(t1h,'Value',num2str(t1));
set(t2h,'Value',num2str(t2));
set(dh,'Value',num2str(dd));

eval(t1h.Callback);
eval(t2h.Callback);
eval(dh.Callback);

Step1_Handle = getSimulinkBlockHandle('Robotics_HW8/Step 1');
Step2_Handle = getSimulinkBlockHandle('Robotics_HW8/Step 2');
Step3_Handle = getSimulinkBlockHandle('Robotics_HW8/Step 3');

set_param(Step1_Handle,'After',num2str(deg2rad(T1)));
set_param(Step2_Handle,'After',num2str(deg2rad(T2)));
set_param(Step3_Handle,'After',num2str(DD/1000));

end
