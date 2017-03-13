% Reads from file
M1=dlmread('Data_1.txt');
M2=dlmread('Data_2.txt');
M4=dlmread('Data_4.txt');
M8=dlmread('Data_8.txt');
M16=dlmread('Data_16.txt');


n=M1(:,1)';
error=M1(:,2)';
figure(1)
loglog(n,error)
hold on
for i =1:4
    p=2^i;
    M=eval(['M', num2str(p)]);
    n=M(:,1)';
    error=M(:,2)';           
    loglog(n,error)
end
title('Loglog plot of error for different n');
ylabel('Error = |S_{n}-\pi|')
xlabel('n')
legend('no threads', '2 threads','4 threads','8 threads','16 threads')


n=M1(:,1)';
time=M1(:,3)';
figure(2)
loglog(n,time)
hold on
for i =1:4
    p=2^i;
    M=eval(['M', num2str(p)]);
    n=M(:,1)';      
    time=M(:,3)';
    loglog(n,time)
end
title('Loglog plot of run time for different n');
ylabel('time in s')
xlabel('n')
legend('no threads', '2 threads','4 threads','8 threads','16 threads')

