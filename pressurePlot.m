clear all
close all
clc

x=dlmread('output/x.dat');
y=dlmread('output/y.dat');
npi=length(x);
npj=length(y);

tend=3600;
dt=0.5;

printTimes=10;
print_dt=printTimes*dt;
printSteps=ceil(tend/print_dt);

fileloc = 'output/p/p_     .00.dat';
count = '    ';

tempMid=zeros(1,printSteps);tempTopRight=tempMid;
tempMean=tempMid;tempMidTop=tempMid;



time=print_dt:print_dt:printSteps*print_dt;
% figure(1)

for n=1:printSteps
    
    fileTime=num2str(n*print_dt);
    
    if length(fileTime)==1
        count(4)=fileTime;
    elseif length(fileTime)==2
        count(3:4)=fileTime;
    elseif length(fileTime)==3
        count(2:4)=fileTime;
    elseif length(fileTime)==4
        count(1:4)=fileTime;
    end
    fileloc(13:16)=count;
    p=dlmread(fileloc);
    
    
    
%     drawnow
%     surf(x(2:npi-1),y(2:npj-1),p(2:npi-1,2:npj-1)')
%     axis([x(2) x(npi-1) y(2) y(npj-1) 0 0.005])
%     
%     colorbar
%     caxis([0 0.005]);
%     
%     F(n)=getframe(gcf);
end

figure('rend','painters','pos',[100 100 900 600])
surf(x(2:npi-1),y(2:npj-1),p(2:npi-1,2:npj-1)')
title(sprintf('t=%g s, n=%g',tend,npi*npj))
axis([x(2) x(npi-1) y(2) y(npj-1) 0 0.005])
xlabel('Width [m]')
ylabel('Height [m]')
zlabel('Pressure [Pa]')
caxis([0 0.005]);
c=colorbar;
c.Label.String = 'Pressure [Pa]';