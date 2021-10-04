rng(0);


num_trials=10;  %no. of runs for each (n,k,matrix) pair.

q = 101;        %prime q
M2 = [9000 13800];  %Values of m for different values of k, in the case of Gaussian Random Matrix


n = q^2;
num_k=2;    %no. of different sparsities to be tested
K = zeros(num_k,1);     
K(1) = 17;          %k=17 and 30 to be tested.
K(2) = 30;

Error = zeros(num_k,2);     %reconstruction errors

T_ldpc = zeros(num_k);      %time taken for reconstruction with ldpc matrix
T_gauss = zeros(num_k);     %time taken for reconstruction with random gaussian matrix
    
for j=1:num_k

    l = K(j)+1;             %calculating m1 (m for ldpc matrix)
    m1 = l*q;
    m2 = M2(j);      
    
    A1=LDPC(m1,n);          %Generating the LDPC matrix, with the given dimensions.

    A2=1/sqrt(m2)*randn(m2,n);  %Random Gaussian Matrix

    for ii=1:num_trials

        x = zeros(n,1); %sparse vector 
        x(randperm(n,K(j)))=2*rand(K(j),1)-1; %random bounded vector
        y1 = A1*x;  %measurements
        y2 = A2*x;

        tic;
        rec1_omp = omp(A1,y1,1e-3);    %omp(A,y,eps)
        End = toc;
        T_ldpc(j) = T_ldpc(j) + End;    %calculating total time for all runs.

        err = sum((rec1_omp-x).^2)/sum(x.^2);
        Error(j,1) = Error(j,1) + sqrt(err);    %summed error for all runs.

        tic;
        rec2_omp = omp(A2,y2,1e-4);
        End = toc;
        T_gauss(j) = T_gauss(j) + End;

        err = sum((rec2_omp-x).^2)/sum(x.^2);
        Error(j,2) = Error(j,2) + sqrt(err);
    end


end
    

Error = Error./num_trials;  %avg error
T_ldpc = T_ldpc./num_trials; %avg time taken
T_gauss = T_gauss./num_trials;

fprintf('The Results for LDPC Matrix are-\n For k=17\n');
fprintf('\t avg Error(RMSE)= %d\n \t avg Time(in sec)= %d\n',Error(1,1),T_ldpc(1));
fprintf('For k=30\n');
fprintf('\t avg Error(RMSE)= %d\n \t avg Time(in sec)= %d\n',Error(2,1),T_ldpc(2));

fprintf('The Results for Random Gaussian Matrix are-\n For k=17\n');
fprintf('\t avg Error(RMSE)= %d\n \t avg Time(in sec)= %d\n',Error(1,2),T_gauss(1));
fprintf('For k=30\n');
fprintf('\t avg Error(RMSE)= %d\n \t avg Time(in sec)= %d\n',Error(2,2),T_gauss(2));

