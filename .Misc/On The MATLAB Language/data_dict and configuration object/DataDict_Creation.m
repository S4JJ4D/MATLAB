%% Variable Definition

n = 5; %system Order
m = 3;%system input order
Ap = -1*eye(n); %system dynamics matrix
Am = -3*eye(n); %parameterized system dynamics matrix
Kp = randi(5,n,m); %system dynamics matrix

%Solving Lyapunov Equation
Q = eye(n);
P = lyap(Am.',Q);


%% Initialization

%Create a data dictionary and get DesignDataSection object
if ~exist('Model_2_vector_DD.sldd','file')
    vectorCaseDDict_obj = Simulink.data.dictionary.create('Model_2_vector_DD.sldd');  % Simulink.data.Dictionary class
    
    %Retrieving design data section:    
    dDataSect_obj = getSection(vectorCaseDDict_obj,'Design Data');  % Simulink.data.dictionary.Section class
    %Adding variables to design section:
    addEntry(dDataSect_obj,'m',m);
    addEntry(dDataSect_obj,'n',n);
    addEntry(dDataSect_obj,'Ap',Ap);
    addEntry(dDataSect_obj,'Am',Kp);
    addEntry(dDataSect_obj,'Kp',Kp);
    addEntry(dDataSect_obj,'P',P);
    %Save changes
    saveChanges(vectorCaseDDict_obj);
    
end

%% Modifing Existing Variables In "Design Data" Section:

assignin(dDataSect_obj,'m',m);
assignin(dDataSect_obj,'n',n);
assignin(dDataSect_obj,'Ap',Ap);
assignin(dDataSect_obj,'Am',Am);
assignin(dDataSect_obj,'Kp',Kp);
assignin(dDataSect_obj,'P',P);
%Save changes
saveChanges(vectorCaseDDict_obj);

%% Retrieve entries in data dictionary

Kp_Entry = dDataSect_obj.getEntry('Kp');  % Simulink.data.dictionary.Entry class
Kp_Entry_Val = Kp_Entry.getValue;  % Simulink.Paramter

%% Configure Model Configuration

% The Configuration is a set of parameters that define settings for a model's execution (simulation) and/or deployment (code generation).
% A Configuration can also be a standalone object without attaching to any model. Users may backup/restore/transfer configuration from one model to another.
%create a configuration set;
ConfigSet = Simulink.ConfigSet;

% to list all parameters for configSet object, use either one of the
% following methods:

props = get_param(ConfigSet, 'ObjectParameters');
props = ConfigSet.getprop;
% Note that you can use methods(ConfigSet) to find more methods.
% second call is recommended. (gives the complete list)

%Set ConfigSet properties
set_param(ConfigSet,"name","CustomConfigSet");
set_param(ConfigSet,"Description","Fixed-Step,ode3 Solve");
% frequently used parameters:
% -Fixed Step, -LoadExternalInput, -ExternalInput, -LoadInitialState, -SaveFormat,
% -SaveOutput, -SaveState, -SaveTime, -SignalLogging, -SimUserSources
% -SolverName, -SolverType, -StartTime, -StopTime, -LoadInitialState,
% -InitialState, -SaveFinalState, -SaveOperatingPoint

%Get Solver Component (This is not necessary, in the sense that, you can ...
%directly access all properties of ConfigSet including solver properties...,
%but this practice provides more modular access to properties and enhances
%readability);

% Also note that getComponent is not documented in help. you can find other
% similar methods useing methods(ConfigSet)
ConfigSet_Solver_Component = getComponent(ConfigSet,'Solver');
% other components: 'Data Import/Export', ...

%Configure solver settings:
% You can access properties of solver component thru get_param(ConfigSet_Solver_Component, 'ObjectParamters')
set_param(ConfigSet_Solver_Component,'StopTime','13',...
                          'SolverType','Fixed-step',...
                          'Solver','auto',...
                          'FixedStep','.001');

                      
% Add the newly created configuration object to data dictionary's
% configuartion section:
CSect_obj = getSection(vectorCaseDDict_obj,'Configurations');
addEntry(CSect_obj,'CustomConfigSet',ConfigSet);

%Defining Configuration set reference object
ConfigSetRef = Simulink.ConfigSetRef;

%Assign a configuration set to ConfigSetRef object
ConfigSetRef.Name = 'ConfigRef_1';
ConfigSetRef.SourceName = 'CustomConfigSet';

%Add the configuration set reference to the model
attachConfigSet(gcs, ConfigSetRef);
%Activate configuration set reference 
setActiveConfigSet(gcs, ConfigSetRef.Name);
% And you can also remove the old ConfigSet. Note that the default name for
% a ConfigSet object is 'Configuration'
detachConfigSet(gcs,'Configuration');

