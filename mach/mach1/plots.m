% Reads from file
M2=dlmread('Data_2_1.txt');
M4=dlmread('Data_4_1.txt');
M8=dlmread('Data_8_1.txt');
M16=dlmread('Data_16_1.txt');
M32=dlmread('Data_32_1.txt');
%First column is different n
% n=M(:,1)';
% %Second column is error correspopnding to different n
% error2=M2(:,2)';
% error4=M4(:,2)';
% error8=M8(:,2)';
% error16=M16(:,2)';
% error32=M32(:,2)';
% %Third column is computational time corresponding to different n
% time2=M2(:,3)';
%time4=
% making n, time and the error vec for the different number of processes

n=M2(:,1)';
error=M2(:,2)';
figure(1)
 loglog(n,error)
hold on
for i =1:5
    p=2^i;
    M=eval(['M', num2str(p)]);
    n=M(:,1)';
    error=M(:,2)';           
    loglog(n,error)
end
title('Loglog plot of error for different n');
ylabel('Error = |S_{n}-\pi|')
xlabel('n')
legend('2 processes', '4 processes','8 processes','16 processes','32 processes')


n=M2(:,1)';
time=M2(:,3)';
figure(2)
loglog(n,time)
hold on
for i =1:5
    p=2^i;
    M=eval(['M', num2str(p)]);
    n=M(:,1)';
    error=M(:,2)';           
    time=M(:,3)';
    loglog(n,time)
end
title('Loglog plot of run time for different n');
ylabel('time in s')
xlabel('n')
legend('2 processes', '4 processes','8 processes','16 processes','32 processes')

