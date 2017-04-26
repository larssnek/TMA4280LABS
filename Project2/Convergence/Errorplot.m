%% Speed up matrix
% S is a mxn matrix where m is number of different processes (except 1) and n is number of n
n=15-4;
m=5;
S=zeros(m,n);

%% Plot error

%% P=1
M=dlmread('maxerror_p1.txt');
n=M(:,1);
error=M(:,2);
time1=M(:,3);
figure(1)
hold on
loglog(n, time1,'*-') %, n, n.^(-2))
hold off
figure(2)
loglog(n, error,'*-', n, n.^(-2))
title('loglog plot of error with p=1')
ylabel('||u_{approx}-u_{exact}||_{\infty}')
xlabel('n')
legend("Max norm of error", "slope 2")

%% P=2
M=dlmread('maxerror_p2.txt');
%M=dlmread('maxerror_p2_lille.txt');
n=M(:,1);
error=M(:,2);
time=M(:,3);
S(1,:)=time1./time;
figure(1)
hold on
loglog(n, time,'*-') %, n, n.^(-2))
hold off
figure()
loglog(n, error,'*-', n, n.^(-2))
title('loglog plot of error with p=2')
ylabel('||u_{approx}-u_{exact}||_{\infty}')
xlabel('n')
legend("Max norm of error", "slope 2")

%% P=4
M=dlmread('maxerror_p4.txt');
n=M(:,1);
error=M(:,2);
time=M(:,3);
S(2,:)=time1./time;
figure(1)
hold on
loglog(n, time,'*-') %, n, n.^(-2))
hold off
figure()
loglog(n, error,'*-', n, n.^(-2))
title('loglog plot of error with p=4')
ylabel('||u_{approx}-u_{exact}||_{\infty}')
xlabel('n')
legend("Max norm of error", "slope 2")

%% P=8
M=dlmread('maxerror_p8.txt');
%M=dlmread('maxerror_p8_lille.txt');
n=M(:,1);
error=M(:,2);
time=M(:,3);
S(3,:)=time1./time;
figure(1)
hold on
loglog(n, time,'*-') %, n, n.^(-2))
hold off
figure()
loglog(n, error,'*-', n, n.^(-2))
title('loglog plot of error with p=8')
ylabel('||u_{approx}-u_{exact}||_{\infty}')
xlabel('n')
legend("Max norm of error", "slope 2")


%% P=16
M=dlmread('maxerror_p16.txt');
% M=dlmread('maxerror_p16_lille.txt');
n=M(:,1);
error=M(:,2);
time=M(:,3);
S(4,:)=time1./time;
figure(1)
hold on
loglog(n, time,'*-') %, n, n.^(-2))
hold off
figure()
loglog(n, error,'*-', n, n.^(-2))
title('loglog plot of error with p=16')
ylabel('||u_{approx}-u_{exact}||_{\infty}')
xlabel('n')
legend("Max norm of error", "slope 2")



%% P=32
M=dlmread('maxerror_p32.txt');
% M=dlmread('maxerror_p16_lille.txt');
n=M(:,1);
error=M(:,2);
time=M(:,3);
S(5,:)=time1./time;
figure(1)
hold on
loglog(n, time,'*-') %, n, n.^(-2))
hold off
figure()
loglog(n, error,'*-', n, n.^(-2))
title('loglog plot of error with p=32')
ylabel('||u_{approx}-u_{exact}||_{\infty}')
xlabel('n')
legend("Max norm of error", "slope 2")


%% Axes timings
figure(1)
title('loglog plot run time')
ylabel('t')
xlabel('n')
legend('2 procs', '4 procs', '8 procs', '16 procs', '32 procs')

%% Speed up


%% Hybrid model
M=dlmread('maxerror_p32.txt');
% M=dlmread('maxerror_p16_lille.txt');
n=M(:,1);
error=M(:,2);
time=M(:,3);
p=[1,2,4,8,16,32];
t=[32,16,8,4,2,1];


