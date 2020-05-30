%% task 1
m = 10;
n = 6;
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(m,n);
b = ones(m,1);
T1_result = zeros(m*n,n);

Q = eye(m);
R = A;
for i = 1:n
    w = R(i:end,i);
    w(1) = w(1) + sign(w(1)) * norm(w);
    w = w/norm(w);
    R(i:end,i:end) = R(i:end,i:end) - 2*w*(w'*R(i:end,i:end));
    Q(:,i:end) = Q(:,i:end) - 2*(Q(:,i:end)*w)*w';
    T1_result((i-1)*m+1:i*m,:) = R;
end
save CP4_T1.dat T1_result -ascii

R3 = [A,b];
for i = 1:n+1
    w = R3(i:end,i);
    w(1) = w(1) + sign(w(1)) * norm(w);
    w = w/norm(w);
    R3(i:end,i:end) = R3(i:end,i:end) - 2*w*(w'*R3(i:end,i:end));
end
x = Backsub(R3(1:n,:));

T1_min = norm(A*x-b);
save CP4_min.dat T1_min -ascii


%% Task 2
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,6);
b = r.randn(10,1);

An = A'*A;
bn = A'*b;
x0 = An\bn;

save CP4_T2.dat x0 -ascii


%% Task 3
m = 30;
n = 10;
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(m,n);
b = r.randn(m,1);
epsilon = 1e-8;
y = zeros(n,1);
c = A'*b;
T3_result = [];
for i = 1:100
    res = A'*(A*y)-c; 
    v = A*res;
    step = (res'*res)/(v'*v)*res;
    y = y - step;
    T3_result = [T3_result, norm(step, Inf)];
    if max(abs(step)) < epsilon
        break;
    end
end

save CP4_T3.dat T3_result -ascii


%% Task 4
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(6,10);
b = r.randn(6,1);
A_new = [(1e-8) * eye(10,10); A];
b_new = [zeros(10,1);b];

R = [A_new,b_new];
for i = 1:10+1
    w = R(i:end,i);
    w(1) = w(1) + sign(w(1)) * norm(w);
    w = w/norm(w);
    R(i:end,i:end) = R(i:end,i:end) - 2*w*(w'*R(i:end,i:end));
end
x = Backsub(R(1:n,:));

save CP4_T4.dat x -ascii





