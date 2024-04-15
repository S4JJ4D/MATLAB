Base_handle = getSimulinkBlockHandle('Robotics_HW8/Base')
Arm1_handle = getSimulinkBlockHandle('Robotics_HW8/Arm_1')
Arm2_handle = getSimulinkBlockHandle('Robotics_HW8/Arm_2')
Manipulator_handle = getSimulinkBlockHandle('Robotics_HW8/Manipulator')
get_param(Base_handle,'ExtGeomFileName')
get_param('Robotics_HW8','FileName')
modelDir = get_param('Robotics_HW8','FileName')
[parentDir,~,~] = fileparts(modelDir)

Arm1_Path = [parentDir '\Arm 1.STEP']
Arm2_Path = [parentDir '\Arm 2.STEP']
Base_Path =  [parentDir '\Base.STEP']
Manipulator_Path = [parentDir '\Manipulator.STEP']

if exist(Base_Path) && exist(Arm1_Path) && exist(Arm2_Path) &&  exist(Manipulator_Path)

set_param(Base_handle,'ExtGeomFileName',Base_Path)
set_param(Arm1_handle,'ExtGeomFileName',Arm1_Path)
set_param(Arm2_handle,'ExtGeomFileName',Arm2_Path)
set_param(Manipulator_handle,'ExtGeomFileName',Manipulator_Path)
msg_handle = msgbox({'STEP files successfully added to the model.','Postload operation Completed.'},'Data Loaded')

else

[FileName, PathName] = uigetfile('*.STEP','Select the file named "Base.STEP" in order to import to the model');
Base_Path = [PathName FileName]
[FileName, PathName] = uigetfile('*.STEP','Select the file named "Arm 1.STEP" in order to import to the model');
Arm1_Path = [PathName FileName]
[FileName, PathName] = uigetfile('*.STEP','Select the file named "Arm 2.STEP" in order to import to the model');
Arm2_Path = [PathName FileName]
[FileName, PathName] = uigetfile('*.STEP','Select the file named "Manipulator.STEP" in order to import to the model');
Manipulator_Path = [PathName FileName]

set_param(Base_handle,'ExtGeomFileName',Base_Path)
set_param(Arm1_handle,'ExtGeomFileName',Arm1_Path)
set_param(Arm2_handle,'ExtGeomFileName',Arm2_Path)
set_param(Manipulator_handle,'ExtGeomFileName',Manipulator_Path)
msg_handle = msgbox({'STEP files successfully added to the model.','Postload Operation Completed.'},'Data Successfully Loaded')


end