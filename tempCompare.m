clear all
close all
clc

tend=3600;
dt=1;

printTimes=10;
print_dt=printTimes*dt;
printSteps=ceil(tend/print_dt);
time=0:print_dt:printSteps*print_dt;

x=dlmread('output/x.dat');
y=dlmread('output/y.dat');

load('tempMeanBTCS.mat');
tempMeanBTCS=tempMean;
load('tempMeanBOUS.mat'); 
tempMeanBOUS=tempMean;
expfile=dlmread('tempEXPERIMENTAL.dat');
tempEXP=interp1(expfile(:,1),expfile(:,2),time);

figure('rend','painters','pos',[100 100 900 600])
hold on
plot(time,tempMeanBTCS,'-','LineWidth',2)
plot(time,tempMeanBOUS,'--','LineWidth',2)
plot(expfile(:,1),expfile(:,2),'+','LineWidth',2)
title('Mean temperature')
axis([0 3600 20 90 ])
xlabel('Time [s]')
ylabel('Temperature [K]')
grid minor
legend('BTCS','Boussinesq','Experimental')