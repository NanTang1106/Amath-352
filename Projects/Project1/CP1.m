%% Task 1
samp_size = [100, 1000, 10000];
samp_results = [];
r = RandStream('mt19937ar','Seed',1234);

for i = 1:3
    total = 0;
    for j = 1:samp_size(i)
        x = r.rand(1);
        y = r.rand(1);
        if y < x^2 
            total = total + 1;
        end
    end
    output = total / samp_size(i);
    samp_results = [samp_results output];
end

save CP1_T1.dat samp_results -ascii


%% Task 2
n = 6;
m = 7;
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(n, m);
A_RGE = zeros(n * (n - 1), m);

for i = 1:(n-1)
    if A(i,i) == 0
            disp("Method failed: matrix is rank deficient")
            break
    else 
        for j = i+1:n
            l_j = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - l_j * A(i,:);
        end
    end
    A_RGE(1+(i-1)*n:i*n,:) = A;
end

save CP1_T2.dat A_RGE -ascii
            

%% Task 3
U = A;
x = U(:,n+1);
U_BWS = zeros(n, n-1);

if U(n,n) == 0
    disp("Method failed: singular matrix")
else 
    x(n) = U(n,n+1)/U(n,n);
    for i = (n-1):-1:1
        if U(i,i) == 0
            disp("Method failed: singular matrix")
            break
        else
            r_sum = 0;
            for j = i+1:n
                r_sum = r_sum + U(i,j) * x(j);
            end
            x(i) = (U(i, n+1) - r_sum) / U(i,i);
        end
        U_BWS(:,n-i) = x;
    end
end

save CP1_T3.dat U_BWS -ascii
        

%% Task 4
n = 6;
m = 7;
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(n, m);
A_GEWP = zeros(n * (n - 1), m);

for i = 1:(n-1)
    r_first = i;
    for k = i:n
        if A(k,i) ~= 0 && abs(A(k,i)) > abs(A(r_first,i))
            r_first = k;
        end
    end
    if A(r_first,i) == 0
        disp("Method failed: singular matrix")
        break
    end
    r_temp = A(i,:);
    A(i,:) = A(r_first,:);
    A(r_first,:) = r_temp;
    for j = i+1:n
        l_j = A(j,i)/A(i,i);
        A(j,:) = A(j,:) - l_j * A(i,:);
    end
    A_GEWP(1+(i-1)*n:i*n,:) = A;
end

save CP1_T4.dat A_GEWP -ascii



