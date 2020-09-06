function[dS1,uc_dS1,dS_t1,dS2,uc_dS2,dS_t2] = calmean(S,uc,t,first,last)

t1 = find(t == first);
t2 = find(t == last);

if first < datenum(2017,6,15) && last > datenum(2018,6,15)
    t3 = find(t == datenum(2017,6,15));
    t4 = find(t == datenum(2018,6,15));
    S_t = S(t1:t2);
    uc_t = uc(t1:t2);
    n1 = t3-t1+1;
    n2 = t2-t4+1;
    n = t2 - t1 + 1;
    A = zeros(n-4,n);
    B = zeros(n-4,n);
    for i = 1 : n1-2
        A(i,i) = -1/2;
        A(i,i+2) = 1/2;
        B(i,i) = 1/4;
        B(i,i+2) = 1/4;
    end
    for i = n1-1 : n-4
        A(i,i+2) = -1/2;
        A(i,i+4) = 1/2;
        B(i,i+2) = 1/4;
        B(i,i+4) = 1/4;
    end
    dS = A * S_t;
    dS1 = dS(1:n1-2);
    dS2 = dS(n1-1:n-4);
    dS_t1 = t(t1+1:t3-1);
    dS_t2 = t(t4+1:t2-1);
    dS_t = [dS_t1;dS_t2];
else
    S_t = S(t1:t2);
    uc_t = uc(t1:t2);
    n = t2 - t1 +1;
    A = zeros(n-2,n);
    B = zeros(n-2,n);
    for i = 1:n-2  
        A(i,i) = -1/2;
        A(i,i+2) = 1/2;
        B(i,i) = 1/4;
        B(i,i+2) = 1/4;
    end
    dS = A * S_t;
    dS_t = t(t1+1:t2-1);
end
uc_dS = sqrt(B * uc_t);
uc_dS1 = uc_dS(1:n1-2);
uc_dS2 = uc_dS(n1-1:n-4);
end