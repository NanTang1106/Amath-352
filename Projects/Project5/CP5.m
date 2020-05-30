%% task 1
gamma = 1;
b = [5; 2.5; 1.5; 5];
A = [1 0.1 1 0.1 1 1; 0.5 0 0.5 0.4 0.5 1;
    0.1 0.1 0.1 0.1 0.1 0.1; 0 1 0 0 1 0.1];
An = [sqrt(gamma) * eye(6); A];
bn = [zeros(6,1); b];
dx = (An'* An)\(An' * bn);
save CP5_T1.dat dx -ascii


%% task 2
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,10);
trueEig = eigs(A);
A = A+A';
i = 0;
err = norm(A - diag(A).*eye(10),1);
while err > 1e-3
    i = i+1;
    [Q,R] = qr(A);
    A = R*Q;
    err = norm(A - diag(A).*eye(10),1);
end 
estEig = diag(A);
save CP5_T2_i.dat i -ascii
save CP5_T2_lam.dat estEig -ascii


%% task 3
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,10)+1j*r.randn(10,10);
[V,D] = eig(A);
diagD = diag(D);
for i = 1:10
    localMin = min(real(diagD(i:10)));
    indexMin = find(real(diagD) == localMin);
    tempV = V(:,indexMin);
    tempRD = diagD(indexMin);
    tempD = D(indexMin,:);
    V(:,indexMin) = V(:,i);
    V(:,i) = tempV;
    diagD(indexMin) = diagD(i);
    diagD(i) = tempRD;
end
D = diagD.* eye(10);
vr = real(V);
vi = imag(V);
dr = real(D);
di = imag(D);

save CP5_T3_Vr.dat vr -ascii
save CP5_T3_Vi.dat vi -ascii
save CP5_T3_Dr.dat dr -ascii
save CP5_T3_Di.dat di -ascii


%% Task 4
%{
url = 'http://faculty.washington.edu/trogdon/352/CP5_M1.mat';
filename = 'CP5_M1.mat';
outfilename = websave(filename,url);
load CP5_M1.mat
url = 'http://faculty.washington.edu/trogdon/352/CP5_M2.mat';
filename = 'CP5_M2.mat';
outfilename = websave(filename,url);
load CP5_M2.mat
url = 'http://faculty.washington.edu/trogdon/352/CP5_M3.mat';
filename = 'CP5_M3.mat';
outfilename = websave(filename,url);
load CP5_M3.mat

A = M2*M3*M1';
%A = rgb2gray(A);
imshow(A)
[U,S,V] = svd(double(A));
%S(1:5,1:5)
k = 1; % A = U S V^T
W = V';
[m,n] = size(A);
B =  U(:,1:k)*S(1:k,1:k)*W(1:k,:);
B = uint8(B);
imshow(B)
fprintf('Compression ratio = %.5f', (m*k+k+n*k)/(m*n))

k = 20; % A = U S V^T
W = V';
[m,n] = size(A);
B =  U(:,1:k)*S(1:k,1:k)*W(1:k,:);
B = uint8(B);
imshow(B)
fprintf('Compression ratio = %.5f', (m*k+k+n*k)/(m*n))
%}

t4 = 17;
save CP5_T4.dat t4 -ascii




