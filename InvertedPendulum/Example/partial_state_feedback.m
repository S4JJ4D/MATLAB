%% 
% s-domain SISO system is converted into SS representation with "physically
% meaningfull" states. We then use partial state feedback to stabilize the
% dyanmics around the origin. 



% +------------------------------------------------------------------+-------------------+
% |                       State-Space Control                        | Classical Control |
% +------------------------------------------------------------------+-------------------+
% | full-state feedback                                              | PD  Control       |
% | Partial state feedback (1st state)                               | P   Control       |
% | Partial state feedback (augmented state + 1st state)             | PI  Control       |
% | full state feedback (augmented state + 1st state + 2nd state)    | PID Control       |
% +------------------------------------------------------------------+-------------------+



% mass-damper-spring system
sys_tf = tf(1, [1 1 1]);
sys_ss = ss([0 1;-1 -1], [0;1], [1 0], 0);
% sys_ss.StateName = {'p', 'v'};

% attempting to place eigenvalues of the system on the LHS of the complex
% plane by partial state feedback:

syms k1 k2 real;
% using only the first state, x1 : pos
cl_sys = sys_ss.A + sys_ss.B * -[k1 0];
vpa(charpoly(cl_sys), 3)
% results in [1.0, 1.0, k1 + 1.0]. It is evident that any k1 > 0 is
% acceptable: both eigenvalues of the close loop system would lie in LHS.
% This is precisely the proportial controller for the original SISO system.

% example:
k1 = 5;
roots([1 1 k1+1])
bh = getSimulinkBlockHandle('classical_modern_interplay/Example_1');
set(bh, 'Commented', 'off'); % uncomment the block
simout = sim('classical_modern_interplay.slx');
time = simout.tout;
figure;
plot(time, simout.logsout.getElement('y_ss').Values.Data, ...
    'b-', 'DisplayName', 'y_{ss}');
hold on;
plot(time, simout.logsout.getElement('y').Values.Data, ...
    'ro', 'DisplayName', 'y');
legend;
set(bh, 'Commented', 'on'); % comment the block

%% Augment integral of a state and use partial state feedback
% | Partial state feedback (augmented state + 1st state) | PI Control        |

% augment integral of the first state (position) to design  PI control

A2 = [[0 1 0];[[0;0],[0 1;-1 -1]]];
B2 = [0;[0;1]];
C2 = [0,[1 0]];
D2 = 0;
sys2_ss = ss(A2,B2,C2,D2);
sys2_ss.StateName = {'z', 'x1', 'x2'};
rank(ctrb(sys2_ss))

syms k1 k2 k3 real;
% using only the first state, x1 : pos
cl_sys = sys2_ss.A + sys2_ss.B * -[k1 k2 0];
vpa(charpoly(cl_sys), 3)
% results in [1.0, 1.0, k2 + 1.0, k1]. Routh-Hurwitz criterion requires
% k2+1 > k1 > 0.

% example:
k1 = .46;
k2 = .68;
roots([1.0, 1.0, k2 + 1.0, k1])
bh = getSimulinkBlockHandle('classical_modern_interplay/Example_2');
set(bh, 'Commented', 'off'); % uncomment the block
simout = sim('classical_modern_interplay.slx');
time = simout.tout;
figure;
plot(time, simout.logsout.getElement('y_ss').Values.Data, ...
    'b-', 'DisplayName', 'y_{ss}');
hold on;
plot(time, simout.logsout.getElement('y').Values.Data, ...
    'ro', 'DisplayName', 'y');
legend;
set(bh, 'Commented', 'on'); % comment the block

%% % | full state feedback (augmented state + 1st state + 2nd state)    | PID Control       |

% i, p, d
K = place(sys2_ss.A, sys2_ss.B, [-2 -3 -4]); 
[k1, k2, k3] = deal(K(1), K(2), K(3));

bh = getSimulinkBlockHandle('classical_modern_interplay/Example_3');
set(bh, 'Commented', 'off'); % uncomment the block
simout = sim('classical_modern_interplay.slx');
time = simout.tout;
figure;
plot(time, simout.logsout.getElement('y_ss').Values.Data, ...
    'b-', 'DisplayName', 'y_{ss}');
hold on;
plot(time, simout.logsout.getElement('y').Values.Data, ...
    'ro', 'DisplayName', 'y');
legend;
set(bh, 'Commented', 'on'); % comment the block

