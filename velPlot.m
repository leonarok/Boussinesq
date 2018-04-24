clear all
close all
clc

x=dlmread('output/x.dat');
y=dlmread('output/y.dat');
npi=length(x);
npj=length(y);

tend=2000;
dt=5;

printTimes=10;
print_dt=printTimes*dt;
printSteps=ceil(tend/print_dt);

fileloc1 = 'output/u/u_     .00.dat';
fileloc2 = 'output/v/v_     .00.dat';
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
    fileloc1(13:16)=count;
    fileloc2(13:16)=count;
    u=dlmread(fileloc1);
    v=dlmread(fileloc2)';
    velTot=sqrt(u.^2+v.^2);
    
%     tempMid(n)=T(npi/2,npj/2);
%     tempMidTop(n)=T(npi/2,ceil(npj*3/4));
%     tempMean(n)=mean(mean(T(2:npi-1,2:npj-1)));
%     tempTopRight(n)=T(ceil(npi*3/4),ceil(npj*3/4));
    
    
    drawnow
%     quiver(x(2:npi-1),y(2:npj-1),u(2:npi-1,2:npj-1)',v(2:npi-1,2:npj-1)')

    %streamslice(x,y,u',v')
    %axis([x(2) x(npi-1) y(2) y(npj-1) 0 0.005])
    %colorbar
    %F(n)=getframe(gcf);
end

figure('rend','painters','pos',[100 100 900 600])
hold on
% pcolor(x,y,velTot')
% colormap(summer)
shading interp
h=quiver(x,y,u',v',1.2,'filled');
set(h,'Color','b');
set(h,'LineWidth',0.5)
set(h,'MaxHeadSize',2)
set(h,'AutoScaleFactor',2)
title(sprintf('t=%g s, n=%g',tend,npi*npj))
axis([x(1) x(npi) y(1) y(npj)])
xlabel('Width [m]')
ylabel('Height [m]')
% zlabel('Pressure [Pa]')
% caxis([0 0.00006]);
% c=colorbar;
% c.Label.String = 'Pressure [Pa]';