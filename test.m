%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROBLEM 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note This approach doesn't account for Rotational Inertia'

syms I  m1 m2 M g F1 F2 F3 F
syms q1(t) q2(t) q3(t) q1_d(t) q1_dd(t) q2_d(t) q2_dd(t) q3_d(t) q3_dd(t)
syms a1c a2c a3c a1 a2 a3 
syms l1 l2 ak
%%%%%%%%%%%%%%%%%%%%%%
% Symbol definitions
%%%%%%%%%%%%%%%%%%%%%%
% q(t) - Represents generalized coordinates
% F - Represents generalized force or torque
% q_d(t) - Represents derivative of q wrt to t
% q_dd(t)- Represents second derivative of q wrt to t
% r - Position of COM wrt to inertial(world) frame
% a's - robot link parameters
% I - Rotational Inertia from links

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Lagrangian 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1 = [q2-l1*sin(q1)
      l1*cos(q1)];

r2 = q2;

r3 = 0;
  
% Getting velocities
v1 = diff(r1);  
v2 = diff(r2);


% Kinetic Energy
T1 = 0.5 * transpose(v1) * v1*m1;
T2 = 0.5 * transpose(v2) * v2*M;


T = T1 + T2 ;
T=simplify(T);

% Potential Energy
V1 = m1*g*l1*cos(q1);
v2 = 0

V = V1 + v2
V=simplify(V);

% Lagrangian
L = T-V;
L = simplify(L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Taking derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Derivative of L wrt to q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=subs(L,diff(q1(t),t),q1_d);
L=subs(L,diff(q2(t),t),q2_d);


dL_dq1 = functionalDerivative(L,q1);
dL_dq2 = functionalDerivative(L,q2);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Derivative of L wrt q_dot
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replacing derivatives with symbols
% Matlab can't take functional derivative with
% respect to a derivative 
dL_dq1_dot= functionalDerivative(L,q1_d);
dL_dq2_dot= functionalDerivative(L,q2_d);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Derivative of dl/dq_dot wrt to t
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Taking derivative of dL_dq with respect to t
dt_dL_dq1_dot = diff(dL_dq1_dot,t)
dt_dL_dq2_dot = diff(dL_dq2_dot,t);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %Equation of Motion
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Replacing derivatives with symbols
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Replacing derivatives with symbols to perform substaction
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q1_d(t),t),q1_dd);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q1(t),t),q1_d);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q2_d(t),t),q2_dd);
dt_dL_dq1_dot=subs(dt_dL_dq1_dot,diff(q2(t),t),q2_d);


dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q1_d(t),t),q1_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q1(t),t),q1_d);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q2_d(t),t),q2_dd);
dt_dL_dq2_dot=subs(dt_dL_dq2_dot,diff(q2(t),t),q2_d);




 %Replacing derivatives with symbols to perform substaction





 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %        Finding Equations
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqn1 = simplify(dt_dL_dq1_dot - dL_dq1) == 0
eqn2 = simplify(dt_dL_dq2_dot - dL_dq2) == F 
