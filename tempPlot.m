clear all
close all
clc

npi=42;
npj=82;

tend=600;
dt=0.5;

printTimes=10;
print_dt=printTimes*dt;
printSteps=ceil(tend/print_dt);

fileloc = 'output/temp/temp_     .00.dat';
count = '    ';

tempMid=zeros(1,printSteps);tempTopRight=tempMid;
tempMean=tempMid;tempMidTop=tempMid;

x=dlmread('output/x.dat');
y=dlmread('output/y.dat');


time=print_dt:print_dt:printSteps*print_dt;
figure(1)


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
    
    tempMid(n)=T(npi/2,npj/2);
    tempMidTop(n)=T(npi/2,ceil(npj*3/4));
    tempMean(n)=mean(mean(T(2:npi-1,2:npj-1)));
    tempTopRight(n)=T(ceil(npi*3/4),ceil(npj*3/4));
    
    
    drawnow
    surf(x(2:npi-1),y(2:npj-1),T(2:npi-1,2:npj-1)')
    axis([x(2) x(npi-1) y(2) y(npj-1) 20 83])
    caxis([20 83]);
    colorbar
    shading interp
    %F(n)=getframe(gcf);
end

figure(2)
hold on
plot(time,tempMean,'.-')
plot(time,tempMid,'o-')
plot(time,tempMidTop,'x-')
plot(time,tempTopRight)
hold off
legend('Mean','Middle','Near top','Near top right')

figure('rend','painters','pos',[100 100 900 600])

surf(x(2:npi-1),y(2:npj-1),T(2:npi-1,2:npj-1)')
title(sprintf('Temperature profile at t=%g s, computed with %g grid points',tend,npi*npj))
axis([x(2) x(npi-1) y(2) y(npj-1) 20 83])
xlabel('Width [m]')
ylabel('Height [m]')
zlabel('Temperature [K]')
caxis([20 83]);
c=colorbar;
c.Label.String = 'Temperature [K]';
shading interp
