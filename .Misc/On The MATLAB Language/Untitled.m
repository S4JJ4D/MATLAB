syms x

desired_st_handle_1 = getSimulinkBlockHandle('Robotics_HW8/DC 1/Compensator/Desired Settling Time');
desired_overshoot_handle_1 = getSimulinkBlockHandle('Robotics_HW8/DC 1/Compensator/Desired Maximum Overshoot');
desired_st_handle_2 = getSimulinkBlockHandle('Robotics_HW8/DC 2/Compensator/Desired Settling Time');
desired_overshoot_handle_2 = getSimulinkBlockHandle('Robotics_HW8/DC 2/Compensator/Desired Maximum Overshoot');
desired_st_handle_3 = getSimulinkBlockHandle('Robotics_HW8/DC 3/Compensator/Desired Settling Time');
desired_overshoot_handle_3 = getSimulinkBlockHandle('Robotics_HW8/DC 3/Compensator/Desired Maximum Overshoot');

st1 = get_param(desired_st_handle_1,'Value');
Mp1 = get_param(desired_overshoot_handle_1,'Value');
st2 = get_param(desired_st_handle_2,'Value');
Mp2 = get_param(desired_overshoot_handle_2,'Value');
st3 = get_param(desired_st_handle_3,'Value');
Mp3 = get_param(desired_overshoot_handle_3,'Value');

modelWorkspaceHandle = get_param('Robotics_HW8','modelworkspace');

Jg = getVariable(modelWorkspaceHandle,'Jg');
Jp = getVariable(modelWorkspaceHandle,'Jp');
Jm = getVariable(modelWorkspaceHandle,'Jm');
J_L1 = getVariable(modelWorkspaceHandle,'J_L1');
J_L2 = getVariable(modelWorkspaceHandle,'J_L2');
Bm = getVariable(modelWorkspaceHandle,'Bm');
Kb = getVariable(modelWorkspaceHandle,'Kb');
Km = getVariable(modelWorkspaceHandle,'Km');
r = getVariable(modelWorkspaceHandle,'r');
R = getVariable(modelWorkspaceHandle,'R');
Ra = getVariable(modelWorkspaceHandle,'Ra');

JJ = Jm + Jg + J_L1/r^2;
BB = Bm + Kb*Km/Ra;

w1 = double(solve(st1 == (8*JJ)/(BB+(Km/Ra)*x))); %x=Kd
v1 = double(vpasolve(Mp1 ==...
   exp((-pi*((BB+w1*Km/Ra)/(2*JJ*(sqrt(x*Km/(JJ*Ra))))))/sqrt(1-((BB+w1*Km/Ra)/(2*JJ*(sqrt(x*Km/(JJ*Ra)))))^2))));
w2 = double(solve(st2 == (8*JJ)/(BB+(Km/Ra)*x))); %x=Kd
v2 = double(vpasolve(Mp2 ==...
   exp((-pi*((BB+w2*Km/Ra)/(2*JJ*(sqrt(x*Km/(JJ*Ra))))))/sqrt(1-((BB+w2*Km/Ra)/(2*JJ*(sqrt(x*Km/(JJ*Ra)))))^2))));
w3 = double(solve(st3 == (8*JJ)/(BB+(Km/Ra)*x))); %x=Kd
v3 = double(vpasolve(Mp3 ==...
   exp((-pi*((BB+w3*Km/Ra)/(2*JJ*(sqrt(x*Km/(JJ*Ra))))))/sqrt(1-((BB+w3*Km/Ra)/(2*JJ*(sqrt(x*Km/(JJ*Ra)))))^2))));

Kd_Gain_Handle_1 = getSimulinkBlockHandle('Robotics_HW8/DC 1/Compensator/Kd_Gain');
Kp_Gain_Handle_1 = getSimulinkBlockHandle('Robotics_HW8/DC 1/Compensator/Kp_Gain');
Kd_Gain_Handle_2 = getSimulinkBlockHandle('Robotics_HW8/DC 2/Compensator/Kd_Gain');
Kp_Gain_Handle_2 = getSimulinkBlockHandle('Robotics_HW8/DC 2/Compensator/Kp_Gain');
Kd_Gain_Handle_3 = getSimulinkBlockHandle('Robotics_HW8/DC 3/Compensator/Kd_Gain');
Kp_Gain_Handle_3 = getSimulinkBlockHandle('Robotics_HW8/DC 3/Compensator/Kp_Gain');

set_param(Kd_Gain_Handle_1,'Gain',num2str(w1));
set_param(Kp_Gain_Handle_1,'Gain',num2str(v1));
set_param(Kd_Gain_Handle_2,'Gain',num2str(w2));
set_param(Kp_Gain_Handle_2,'Gain',num2str(v2));
set_param(Kd_Gain_Handle_3,'Gain',num2str(w3));
set_param(Kp_Gain_Handle_3,'Gain',num2str(v3));

if getVariable(modelWorkspaceHandle,'isActivated')

xi = getVariable(modelWorkspaceHandle,'X_i');
yi = getVariable(modelWorkspaceHandle,'Y_i');
zi = getVariable(modelWorkspaceHandle,'Z_i');
xf = getVariable(modelWorkspaceHandle,'X_f');
yf = getVariable(modelWorkspaceHandle,'Y_f');
zf = getVariable(modelWorkspaceHandle,'Z_f');

[t1,t2,dd] = Inv_Kin_Scara(xi,yi,zi);
[T1,T2,DD] = Inv_Kin_Scara(xf,yf,zf);

modelWorkspaceHandle.clear('theta_1','theta_2','d');

modelWorkspaceHandle.assignin('theta_1',t1);
modelWorkspaceHandle.assignin('theta_2',t2);
modelWorkspaceHandle.assignin('d',dd);

Step1_Handle = getSimulinkBlockHandle('Robotics_HW8/Step 1');
Step2_Handle = getSimulinkBlockHandle('Robotics_HW8/Step 2');
Step3_Handle = getSimulinkBlockHandle('Robotics_HW8/Step 3');

set_param(Step1_Handle,'After',num2str(deg2rad(T1)));
set_param(Step2_Handle,'After',num2str(deg2rad(T2)));
set_param(Step3_Handle,'After',num2str(DD/1000));

end
