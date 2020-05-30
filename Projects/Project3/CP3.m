%% Task 1
T1_result = zeros(4,4);
n = 4;
A = [1.1,.2,-.2,.5;
     .2,.9,.5,.3;
     .1,0.,1.,.4;
     .1,.1,.1,1.2];
b = [1;0;1;0];

% Jacobi's
L = tril(A, -1);
U = triu(A, +1);
LpU = L + U;
D = diag(A);
c = b./D;
y = zeros(n,1);
for i = 1:20
    y = - (LpU * y)./D + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-2
        T1_result(1,1) = i;
        T1_result(2,1) = max_err;
        break
    end
end

y = zeros(n,1);
for i = 1:20
    y = - (LpU * y)./D + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-4
        T1_result(1,2) = i;
        T1_result(2,2) = max_err;
        break
    end
end

y = zeros(n,1);
for i = 1:50
    y = - (LpU * y)./D + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-6
        T1_result(1,3) = i;
        T1_result(2,3) = max_err;
        break
    end
end

y = zeros(n,1);
for i = 1:50
    y = - (LpU * y)./D + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-8
        T1_result(1,4) = i;
        T1_result(2,4) = max_err;
        break
    end
end

% Gauss-Seidel's
LpD = tril(A);
c = (LpD)^(-1)*b;

y = zeros(n,1);
for i = 1:20
    y = - (LpD)^(-1)*U*y + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-2
        T1_result(3,1) = i;
        T1_result(4,1) = max_err;
        break
    end
end

y = zeros(n,1);
for i = 1:20
    y = - (LpD)^(-1)*U*y + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-4
        T1_result(3,2) = i;
        T1_result(4,2) = max_err;
        break
    end
end

y = zeros(n,1);
for i = 1:20
    y = - (LpD)^(-1)*U*y + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-6
        T1_result(3,3) = i;
        T1_result(4,3) = max_err;
        break
    end
end

y = zeros(n,1);
for i = 1:20
    y = - (LpD)^(-1)*U*y + c;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-8
        T1_result(3,4) = i;
        T1_result(4,4) = max_err;
        break
    end
end

save CP3_T1.dat T1_result -ascii


%% Task 2
T2_result = zeros(2,2);
x0 = [0.9;0.09;0.01];
p = 0;
M = [1-p-1/200,0,1/10000;
    1/200,999/1000,0;
    p,1/1000,9999/10000];
x = x0;
for i = 1:1000
    x = M * x;
    if x(2) >= 0.5
        T2_result(1,1) = i;
        break
    end
end
x_old = x0;
for i = 1:100000
    x_new = M*x_old;
    max_err = max(abs(x_new - x_old));
    if max_err < 1.0e-16
        T2_result(2,1) = x_new(2);
        break
    end
    x_old = x_new;
end

p = 2/1000;
M = [1-p-1/200,0,1/10000;
    1/200,999/1000,0;
    p,1/1000,9999/10000];
x = x0;
for i = 1:1000
    x = M * x;
    if x(2) >= 0.5
        T2_result(1,2) = i;
        break
    end
end
x_old = x0;
for i = 1:100000
    x_new = M*x_old;
    max_err = max(abs(x_new - x_old));
    if max_err < 1.0e-16
        T2_result(2,2) = x_new(2);
        break
    end
    x_old = x_new;
end

save CP3_T2.dat T2_result -ascii


%% Task 3
k = 40;
Mk = zeros(k+2, k+1);
Mk(1,1) = 1;
Mk(2,1) = 1;
for i = 1:20
    temp_val = (1/(2*i+1));
    temp_val2 = 1/(2*i);
    Mk(1, 2*i + 1) = temp_val;
    Mk(2*(i+1), 2*i+1) = temp_val;
    Mk(1,2*i) = -temp_val2;
    Mk(2*i+1, 2*i) = temp_val2;
end
save CP3_T3_1.dat Mk -ascii

a_vec = zeros(k+1,1);
for i = 0:20
    a_vec(2*i+1,1) = (-1)^i * 1/factorial(i);
end
integ = sum(Mk * a_vec);
save CP3_T3_2.dat integ -ascii


