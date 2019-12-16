function c = compute(m1,m2,M,l1,l2,g)
    syms t
    A = [0 1 0 0 0 0;0 0 m1*g/M 0 m2*g/M 0;0 0 0 1 0 0;0 0 g*m1*((1/l1)+M) 0 g*m2/(l1*M) 0;0 0 0 0 0 1;0 0 g*m1/(l2*M) 0 g*m2*((1/l2)+M) 0]
    B = [0;1/(M);0;1/(M*l1);0;1/(M*l2)]
    eig(A);
    %eigen = [-500 -200 -800 -400 -100 -600];
    %K = place(A,B,eigen)
    Q = 0.0000000000001*eye(6);
    R=0.000000001;
    %expm(A'*t)*Q*expm(A*t)%int(simplify(expm(A'*t)*Q*expm(A*t))
    [K,P,E] = lqr(A,B,Q,R)
    eigen = eig(A-B*K)
    %[K,P,E] = lqr(,B,Q,R)
    c = [B A*B A^2*B A^3*B A^4*B A^5*B];
end