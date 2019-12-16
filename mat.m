function [Af,Bf] = mat(a,b,c,d,e,f)
    syms m1 m2 g M l1 l2 F t1 t2 t1_dot t2_dot x x_dot
    A = (F - m1*g*sin(2*t1)/2 - m2*g*sin(2*t2)/2 - m1*l1*sin(t1)*(t1_dot)^2 - m2*l2*sin(t2)*(t2_dot)^2)/(M + m1*(sin(t1)^2) + m2*(sin(t2)^2));
    B = (1/l1)*(cos(t1)*(F - m1*g*sin(2*t1)/2 - m2*g*sin(2*t2)/2 - m1*l1*sin(t1)*(t1_dot)^2 - m2*l2*sin(t2)*(t2_dot)^2)/(M + m1*(sin(t1)^2) + m2*(sin(t2)^2)) - g*sin(t1));
    C = (1/l2)*(cos(t2)*(F - m1*g*sin(2*t1)/2 - m2*g*sin(2*t2)/2  - m1*l1*sin(t1)*(t1_dot)^2 - m2*l2*sin(t2)*(t2_dot)^2)/(M + m1*(sin(t1)^2) + m2*(sin(t2)^2)) - g*sin(t2));
    f1x1 = diff(A,x);
    f1x1_dot = diff(A,x_dot);
    f1t1 = diff(A,t1);
    f1t1_dot = diff(A,t1_dot);
    f1t2 = diff(A,t2);
    f1t2_dot = diff(A,t2_dot);
    f2x1 = diff(B,x);
    f2x1_dot = diff(B,x_dot);
    f2t1 = diff(B,t1);
    f2t1_dot = diff(B,t1_dot);
    f2t2 = diff(B,t2);
    f2t2_dot = diff(B,t2_dot);
    f3x1 = diff(C,x);
    f3x1_dot = diff(C,x_dot);
    f3t1 = diff(C,t1);
    f3t1_dot = diff(B,t1_dot);
    f3t2 = diff(C,t2);
    f3t2_dot = diff(C,t2_dot);
    
    A11 = subs(f1x1, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A12 = subs(f1t1, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A13 = subs(f1t2, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A14 = subs(f1x1_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A15 = subs(f1t1_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A16 = subs(f1t2_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A21 = subs(f2x1, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A22 = subs(f2t1, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A23 = subs(f2t2, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A24 = subs(f2x1_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A25 = subs(f2t1_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A26 = subs(f2t2_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A31 = subs(f3x1, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A32 = subs(f3t1, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A33 = subs(f3t2, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A34 = subs(f3x1_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A35 = subs(f3t1_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    A36 = subs(f3t2_dot, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});

    f2F = diff(A,F);
    f4F = diff(B,F);
    f6F = diff(C,F);
    B11 = subs(f2F, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    B12 = subs(f4F, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    B13 = subs(f6F, {x,x_dot,t1,t1_dot,t2,t2_dot}, {0,0,0,0,0,0});
    B1 = simplify([0;B11;0;B12;0;B13])
    A1 = simplify([0 1 0 0 0 0; A11 A14 A12 A15 A13 A16;0 0 0 1 0 0; A21 A24 A22 A25 A23 A26;0 0 0 0 0 1;A31 A34 A32 A35 A33 A36])
    c =  simplify(simplify([B1 A1*B1 A1^2*B1 A1^3*B1 A1^4*B1 A1^5*B1]))
    B1
    A1*B1
    A1^2*B1
    A1^3*B1
    A1^4*B1
    A1^5*B1
    det(c)
    %Bf = (subs(B1, {m1,m2,M,l1,l2,g},{a,b,c,d,e,f}))
    Bf = B1
    %Af = (subs(A1, {m1,m2,M,l1,l2,g},{a,b,c,d,e,f}))
    Af = A1
end