clear all
%% task 1
n = 11;
sum_count = 0;
for i = 1:n
    sum_count = sum_count + i;
end

save CP0_T1.dat sum_count -ascii

%% Task 2
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,10);
output = A(7:10, 7:10);

save CP0_T2.dat output -ascii