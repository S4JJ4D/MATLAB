%%
%% Variable Definition
ctrl_param = struct('Kp', 1, 'Kv', 2, 'tau', 1, 'Kd', 4);

%% Initialization

%Create a data dictionary and get DesignDataSection object
Simulink.data.dictionary.closeAll;  % Close all open dictionaries
if ~exist('custom_dict.sldd','file')
    vectorCaseDDict_obj = Simulink.data.dictionary.create('custom_dict.sldd');  % Simulink.data.Dictionary class
    
    %Retrieving design data section:
    dDataSect_obj = getSection(vectorCaseDDict_obj,'Design Data');  % Simulink.data.dictionary.Section class
    %Adding variables to design section:
    addEntry(dDataSect_obj,'ctrl_param',ctrl_param);
    %Save changes
       
    BusObject_Param = Simulink.Parameter;
    BusObject_Param.Value = ctrl_param;
    
    addEntry(dDataSect_obj,'ctrl_param_simparam',BusObject_Param);
    
    saveChanges(vectorCaseDDict_obj);
    
end
