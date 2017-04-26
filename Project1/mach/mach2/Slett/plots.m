% Reads from file
M=dlmread('Data.txt');
%First column is different n
n=M(:,1)';
%Second column is error correspopnding to different n
error=M(:,2)';
%Third column is computational time corresponding to different n
time=M(:,3)';

%Error plot
figure(1);
title('Errorplot for different n');
ylabel('Error = |S_{n}-\pi|')
xlabel('n')
loglog(n, error, '.-')

%Plot of computational time
figure(2);
title('Computational time');
ylabel('Time')
xlabel('n')
loglog(n, time, '.-')