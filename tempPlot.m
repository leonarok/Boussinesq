clear all
close all
clc

npi=42;
npj=82;

tend=3600;
dt=0.5;

printTimes=10;
print_dt=printTimes*dt;
printSteps=ceil(tend/print_dt);

fileloc = 'output/temp/temp_     .00.dat';
count = '    ';

tempMid=83*ones(1,printSteps);tempTopRight=tempMid;
tempMean=tempMid;tempMidTop=tempMid;
tempBotRight=tempMid;

x=dlmread('output/x.dat');
y=dlmread('output/y.dat');


time=0:print_dt:printSteps*print_dt;
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
    fileloc(19:22)=count;
    T=dlmread(fileloc)-273.16;
    
    tempMid(n+1)=T(npi/2,npj/2);
    tempMidTop(n+1)=T(npi/2,ceil(npj*5/6));
    tempMean(n+1)=mean(mean(T(2:npi-1,2:npj-1)));
    tempTopRight(n+1)=T(ceil(npi*5/6),ceil(npj*5/6));
    tempBotRight(n+1)=T(ceil(npi*1/6),ceil(npj*1/6));
    
%     drawnow
%     surf(x(2:npi-1),y(2:npj-1),T(2:npi-1,2:npj-1)')
%     axis([x(2) x(npi-1) y(2) y(npj-1) 20 83])
%     caxis([20 83]);
%     colorbar
%     shading interp

end

figure('rend','painters','pos',[100 100 900 600])
hold on
plot(time,tempMean,'-','LineWidth',2)
plot(time,tempMid,'-','LineWidth',2)
plot(time,tempMidTop,'-','LineWidth',2)
plot(time,tempTopRight,'-','LineWidth',2)
plot(time,tempBotRight,'-','LineWidth',2)
title('Temperature development at selected areas')
axis([0 3600 20 90 ])
xlabel('Time [s]')
ylabel('Temperature [K]')
grid minor
legend('Mean','Middle','Near top','Near top right','Near bottom right')

figure('rend','painters','pos',[100 100 900 600])
surf(x(2:npi-1),y(2:npj-1),T(2:npi-1,2:npj-1)')
title(sprintf('t=%g s, n=%g',tend,npi*npj))
axis([x(2) x(npi-1) y(2) y(npj-1) 20 83])
xlabel('Width [m]')
ylabel('Height [m]')
zlabel('Temperature [K]')
caxis([20 83]);
c=colorbar;
c.Label.String = 'Temperature [K]';
