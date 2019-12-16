function c = compute(m1,m2,M,l1,l2,g)
    syms t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   LQR CONTROLLER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[A,B] = mat(m1,m2,M,l1,l2,g);
    
A = [0,1.0000,0,0,0,0;
         0,0,-1.0000,0,-1.0000,0;
         0,0,0,1.0000,0,0;
         0,0,-0.5500,0,-0.0500,0;
         0,0,0,0,0,1.0000;
         0,0,-0.1000,0,-1.1000,0];
     B = [0;0.001;0;0.00005;0;0.0001];
    
    Eigen = eig(A)
    disp('The eigen values of A ')
    disp('therefore Lyapunovs indirect method is inconclusive for the system')
    tspan = 0:0.1:100;
    Q = 100*eye(6);
    R= 0.01;
    [K,P,E] = lqr(double(A),double(B),Q,R)
    eigen = eig(A-B*K)
    C1 = [1 0 0 0 0 0]              %For output (x(t))
    C2 = [0 0 1 0 0 0;            %For output (theta1(t),theta2(t)) 
          0 0 1 0 1 0]
    C3 = [1 0 0 0 0 0;            %For output (x(t),theta2(t))    
          0 0 0 0 1 0]
    C4 = [1 0 0 0 0 0;             %For output (x(t),theta1(t),theta2(t))
          0 0 1 0 0 0;
          0 0 0 0 1 0]
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %Checking Observability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ob1 = [C1;C1*A;C1*A^(2);C1*A^(3);C1*A^(4);C1*A^(5)]'
    ob2 = [C2;C2*A;C2*A^(2);C2*A^(3);C2*A^(4);C2*A^(5)]'
    ob3 = [C3;C3*A;C3*A^(2);C3*A^(3);C3*A^(4);C3*A^(5)]'
    ob4 = [C4;C4*A;C4*A^(2);C4*A^(3);C4*A^(4);C4*A^(5)]'

    observability1 = rank(ob1)
    disp('Rank is 6. The system is observable for output x(t)')
    observability2 = rank(ob2)
    disp('Rank is 4. The system is not observable for output (t1(t),t2(t))')
    observability3 = rank(ob3)
    disp('The system is observable for output (x(t),t2(t))')
    observability4 = rank(ob4)
    disp('Rank is 6.The system is observable for output (x(t),t1(t),t2(t))')
    s0 = [1; 0; 0; 0; 0; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Obsevability of Linearized System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t,y1] = ode45(@(t,y)(A-B*K)*y,tspan,s0);
    figure;
    hold on
    plot(t,y1(:,1),'r')
    plot(t,y1(:,3),'b')
    plot(t,y1(:,5),'k')
    ylabel('state variables')
    xlabel('time in s')
    title('Response of Linearized system with LQR based control')
    legend('x position of the cart','theta1','theta2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Obsevability of Original Non-Linear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [t1,y2] = ode45(@(t,y)nlinear(y,t,m1,m2,M,l1,l2,g,-K*y),tspan,s0);
    figure;
    hold on
    plot(t1,y2(:,1),'r')
    plot(t1,y2(:,3),'b')
    plot(t1,y2(:,5),'k')
    ylabel('state variables')
    xlabel('time in s')
    title('Response of Non-Linear system with LQR based control')
    legend('x position of the cart','theta1','theta2')
    
    %Introducing Noise and Disturbances in the system
Bd = 0.1*eye(6);             %input disturbance covarianve
Bn = 0.1;
Bn1 = 0;                      %output measurement noise
Bn3 = 0*[0,1;0,1];
Bn4 = 0*[0,0,1;0,0,1;0,0,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Obtaining "best" Luenberger observer for each one of the output vectors using Kalman Bucy Filter (Lqe))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system
[L1,P,E] = lqe(A,Bd,C1,Bd,Bn);
[L3,P,E] = lqe(A,Bd,C3,Bd,Bn*eye(2));
[L4,P,E] = lqe(A,Bd,C4,Bd,Bn*eye(3));
%% Luenberger's Observer using Pole placement Method
%%
Ae=[(A-B*K)];
poles = eig(Ae) 
P = [-2 -5 -6 -7 -8 -9];
L1p = place(A',C1',P)'
L3p = place(A',C3',P)'
C4p = [1,0,0,0,0,0;0,0,1,0,0,0;0,0,0,0,1,0]
L4p = place(A',C4',P)';


% Creating Augmented Matrices for Simulation
uD = randn(6,size(tspan,2));      %input for disturbance
uN = randn(size(tspan));          %input for noise
u = 0*tspan;
u(200:length(tspan)) = 1;      % Step input at t = 10
u1 = [u; Bd*Bd*uD; uN];

uDp = 0*randn(6,size(tspan,2));
uNp = 0*randn(size(tspan));
up = 0*tspan;
up(200:length(tspan)) = 1;      % Step input at t = 10

u1p = [up; Bd*Bd*uDp; uNp];

Be = [B,Bd,zeros(size(B))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Luenberger Observer output when X(t) is the output vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
sysLO1 = ss(A-L1*C1,[B L1],C1,zeros(1,2));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De1 = [0,0,0,0,0,0,0,Bn1];                     %Augmented D matrix

sys1 = ss(A,Be,C1,De1);
[y1,t] = lsim(sys1,u1,tspan);

%Simulating the States of the output variables 
[x1,t] = lsim(sysLO1,[u; y1'],tspan);

figure();
hold on
plot(t,y1(:,1),'g')
plot(t,x1(:,1),'k--')
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C1: output vector x(t) - Linear System')
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout2]=sim('nonlinerLO',tspan,opt);
figure();
hold on
plot(tout2,out1(:,1),'r')
plot(tout2,states1(:,1),'k--')
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C1: output vector x(t) - Nonlinear System')
hold off
%{
%% Luenberger Observer for C1 = [1,0,0,0,0,0] - Pole placement method
%%
sysEstp = ss(A-L1p*C1,[B L1p],C1,zeros(1,2));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De1 = [0,0,0,0,0,0,0,Bn1];                     %Augmented D matrix

sys1 = ss(A,Be,C1,De1);
[y1p,t] = lsim(sys1,u1p,tspan);

%Simulating the States of the output variables 
[x1p,t] = lsim(sysEstp,[up; y1p'],tspan);

figure();
hold on
plot(t,y1p(:,1),'g','Linewidth',2)
plot(t,x1p(:,1),'k--','Linewidth',2)
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C1: output vector x(t) - Linear System')
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout2]=sim('nonlinerLOp',tspan,opt);
figure();
hold on
plot(tout2,out1p(:,1),'r','Linewidth',2)
plot(tout2,states1p(:,1),'k--','Linewidth',2)
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C1: output vector x(t) - Nonlinear System')
hold off
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Luenberger Observer output when (X(t),theta2(t)) is the output vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sysLO3 = ss(A-L3*C3,[B L3],C3,zeros(2,3));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De3 = [zeros(size(C3)),Bn3];                     %Augmented D matrix

sys3 = ss(A,Be,C3,De3);
[y3,t] = lsim(sys3,u1,tspan);

%Simulating the States of the output variables 
[x3,t] = lsim(sysLO3,[u; y3'],tspan);

figure();
hold on
plot(t,y3(:,1),'g')
plot(t,y3(:,2),'b')
plot(t,x3(:,1),'k--')
plot(t,x3(:,2),'r--')
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Estimated Response for C3: output vector (x(t),t2(t))')
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout3]=sim('nonlinearLO3',tspan,opt);
figure();
hold on
plot(tout3,out3(:,1),'r')
plot(tout3,out3(:,2),'g')
plot(t,states3(:,1),'k--')
plot(t,states3(:,2),'r--')
ylabel('x-position of cart')
xlabel('time in s')
legend('Output obtained from noisy system','Estimated output of the system')
title('Estimated Response for C3: output vector (x(t),t2(t)) - Nonlinear System')
hold off
%{
%% Luenberger Observer for C3: For output(x(t),t2(t))-Pole Placement Method
%%
sysLO3p = ss(A-L3p*C3,[B L3p],C3,zeros(2,3));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De3 = [zeros(size(C3)),Bn3];                     %Augmented D matrix

sys3 = ss(A,Be,C3,De3);
[y3p,t] = lsim(sys3,u1p,tspan);

%Simulating the States of the output variables 
[x3p,t] = lsim(sysLO3p,[up; y3p'],tspan);

figure();
hold on
plot(t,y3p(:,1),'g','Linewidth',2)
plot(t,y3p(:,2),'b','Linewidth',2)
plot(t,x3p(:,1),'k--','Linewidth',1)
plot(t,x3p(:,2),'r--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Estimated Response for C3: output vector (x(t),t2(t))-Pole Method' )
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout3p]=sim('nonlinearLO3p',tspan,opt);
figure();
hold on
plot(tout3p,out3p(:,1),'r','Linewidth',2)
plot(tout3p,out3p(:,2),'g','Linewidth',2)
plot(t,states3p(:,1),'k--','Linewidth',1)
plot(t,states3p(:,2),'r--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Estimated Response for C3: output vector (x(t),t2(t))-Pole Method' )
hold off
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Luenberger Observer output when (X(t),theta1(t), theta2(t)) is the output vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sysLO4 = ss(A-L4*C4,[B L4],C4,zeros(3,4));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De4 = [zeros(3,5),Bn4];                     %Augmented D matrix

sys4 = ss(A,Be,C4,De4);
[y4,t] = lsim(sys4,u1,tspan);

%Simulating the States of the output variables 
[x4,t] = lsim(sysLO4,[u;y4'],tspan);

figure();
hold on
plot(t,y4(:,1),'g')
plot(t,y4(:,2),'b')
plot(t,y4(:,3),'c')
plot(t,x4(:,1),'m--')
plot(t,x4(:,2),'r--')
plot(t,x4(:,3),'k--')
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off

opt = simset('solver','ode45','SrcWorkspace','Current');
[tout4]=sim('nonlinearLO4',tspan,opt);

figure();
hold on
plot(tout4,out4(:,1),'g')
plot(tout4,out4(:,2),'b')
plot(tout4,out4(:,3),'c')
plot(tout4,states4(:,1),'m--')
plot(tout4,states4(:,2),'r--')
plot(tout4,states4(:,3),'k--')
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off
%{
%% Luenberger Observer for C4: For output (x(t),t1(t),t2(t)) - Pole Placement Method
%%
sysLO4p = ss(A-(L4p*C4),[B L4p],C4,zeros(3,4));     %State Estimator system

%Obtaining Y values for a system simulated with noise and disturbance. 
De4 = [zeros(3,5),Bn4];                     %Augmented D matrix

sys4 = ss(A,Be,C4,De4);
[y4p,t] = lsim(sys4,u1p,tspan);

%Simulating the States of the output variables 
[x4p,t] = lsim(sysLO4p,[up;y4p'],tspan);

figure();
hold on
plot(t,y4p(:,1),'g','Linewidth',2)
plot(t,y4p(:,2),'b','Linewidth',3)
plot(t,y4p(:,3),'c','Linewidth',2)
plot(t,x4p(:,1),'m--','Linewidth',1)
plot(t,x4p(:,2),'r--','Linewidth',2)
plot(t,x4p(:,3),'k--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off
opt = simset('solver','ode45','SrcWorkspace','Current');
[tout4]=sim('nonlinearLO4p',tspan,opt);
figure();
hold on
plot(tout4,out4p(:,1),'g','Linewidth',2)
plot(tout4,out4p(:,2),'b','Linewidth',3)
plot(tout4,out4p(:,3),'c','Linewidth',2)
plot(tout4,states4p(:,1),'m--','Linewidth',1)
plot(tout4,states4p(:,3),'r--','Linewidth',2)
plot(tout4,states4p(:,5),'k--','Linewidth',1)
ylabel('State Variables ')
xlabel('time in s')
legend('Noisy output x(t)','Noisy output theta1(t)','Noisy output theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Estimated Response for C4: output vector (x(t),t1(t),t2(t))')
hold off
%}
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LQG Controller for smallest Output Vector C1 = [1,0,0,0,0,0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ac = A-L1p*C1;
Bc = [B L1p];
Cc = eye(6);
Dc = 0*[B L1p];

opt = simset('solver','ode45','SrcWorkspace','Current');
sim('nonlinearlqg',tspan,opt);
%% Simulation Results
%%
figure();
hold on
plot(tout,states(:,1),'r')
plot(tout,states(:,3),'k')
plot(tout,states(:,5),'b')
plot(tout,inputlqg(:,1),'g')

title('LQG Nonlinear System')
legend('x-position','theta1','theta2')
hold off


    c = [B A*B A^2*B A^3*B A^4*B A^5*B];
end
