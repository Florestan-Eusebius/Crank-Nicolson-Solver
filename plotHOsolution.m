test=real(readComplexMatrix('HO.txt'));
%%%%
x=linspace(-5,5,501);

box on;
% set(figure,'visible','off');

F(501) = struct('cdata',[],'colormap',[]);
psi0=[test(1,:) test(1,1)];
T=linspace(0,4*pi,501);

for i=1:501
    clf;
    psi=[test(i,:) test(i,1)];
    plot(x,psi,':');
    hold on;
    plot(x,cos(T(i)/2).*psi0,'--');
    ylim([-(1/pi)^0.25 (1/pi)^0.25]);
    xlabel('$x$','Interpreter','Latex');
    ylabel('$\mathrm{Re}\psi(x)$','Interpreter','Latex')
    legend('Numerical result','Theoretical result','location','northeast')
    grid minor;
    F(i)=getframe(gcf);
end

% movie(F);

writerObj =VideoWriter('HO.avi'); % 生成一个avi动画
writerObj.FrameRate=20; % 设置avi动画的参数，设置帧速率
open(writerObj); % 打开avi动画
writeVideo(writerObj,F); % 将保存的动画写入到视频文件中
close(writerObj);

% x=linspace(0,4*pi,501);
% plot(x,test(:,100))
% hold on
% plot(x,test(1,100)*cos(x/2))
% saveas(gca,'test6.png','png')