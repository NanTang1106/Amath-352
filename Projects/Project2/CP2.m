%% Task 1
r = RandStream('mt19937ar', 'Seed', 1234);
n = 6;
m = 6;
A = r.randn(n,m);
P = (1:n)';
M = A;
T1_result = zeros(n*(n-1), (n+1));

for i = 1:n-1
    max_index = i;
    for j = i+1:n
        if abs(M(j,i)) >= abs(M(max_index,i))
            max_index = j;
        end
    end
    P_temp = P(i);
    P(i) = P(max_index);
    P(max_index) = P_temp;
    M_temp = M(i,:);
    M(i,:) = M(max_index,:);
    M(max_index,:) = M_temp;
    if M(i,i) == 0
        disp("Matrix is Rank deficient")
    end
    for j = i+1:n
        a = M(j,i) / M(i,i);
        M(j,i:end) = M(j,i:end) - a .* M(i,i:end);
        M(j,i) = a;
    end
    T1_result(1 + (i-1)*n :i*n,:) = [P, M];
end

save CP2_T1.dat T1_result -ascii


%% Task 2
r = RandStream('mt19937ar', 'Seed', 1234);
n = 6;
m = 6;
A = r.randn(n,m);
b = [1;0;1;0;1;0];

url = 'http://faculty.washington.edu/trogdon/352/CP2_T1.mat';
filename = 'CP2_T1.mat';
outfilename = websave(filename,url);
load CP2_T1.mat

p = O(:,1);
LU = O(:,2:7);
L = tril(LU,-1) + eye(n,n);
U = triu(LU);
z = applyP(p, b);
y = Forsub([L,z]);
x = Backsub([U,y]);
T2_result = [z, y, x];

save CP2_T2.dat T2_result -ascii


%% Task 3
A = [1.1,.2,-.2,.5;
     .2,.9,.5,.3;
     .1,0.,1.,.4;
     .1,.1,.1,1.2];
b = [1;0;1;0];
n = 4;
y = zeros(n,1);
M = eye(n) - A;
T3_result = zeros(3,2);
for i = 1:50
    y = M * y + b;
    max_err = max(abs(A * y - b));   
    if max_err < 1.0e-2
        T3_result(1,1) = i;
        T3_result(1,2) = max_err;
        break
    end
end
y = zeros(n,1);
for i = 1: 50
    y = M * y + b;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-4
        T3_result(2,1) = i;
        T3_result(2,2) = max_err;
        break
    end 
end
y = zeros(n,1);
for i = 1:50
    y = M * y + b;
    max_err = max(abs(A * y - b));
    if max_err < 1.0e-6
        T3_result(3,1) = i;
        T3_result(3,2) = max_err;
        break
    end
end

save CP2_T3.dat T3_result -ascii







