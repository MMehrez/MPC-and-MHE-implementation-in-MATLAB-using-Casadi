function Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;

%--------------------------------------------------------------------------
%-----------------------Simulate robots -----------------------------------
%--------------------------------------------------------------------------
x_r_1 = [];
y_r_1 = [];

r = rob_diam/2;  % obstacle radius
ang=0:0.005:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

figure(500)
% Animate the robot motion
%figure;%('Position',[200 200 1280 720]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

for k = 1:size(xx,2)
    plot([0 12],[1 1],'-g','linewidth',1.2);hold on % plot the reference trajectory
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on % plot exhibited trajectory
    if k < size(xx,2) % plot prediction
        plot(xx1(1:N,1,k),xx1(1:N,2,k),'r--*')
    end
    
    plot(x1,y1,'-sk','MarkerSize',25)% plot robot position
    hold off
    %figure(500)
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    axis([-1 16 -0.5 1.5])
    pause(0.2)
    box on;
    grid on
    %aviobj = addframe(aviobj,gcf);
    drawnow
    % for video generation
    F(k) = getframe(gcf); % to get the current frame
end
close(gcf)
%viobj = close(aviobj)
%video = VideoWriter('exp.avi','Uncompressed AVI');

% video = VideoWriter('exp.avi','Motion JPEG AVI');
% video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
% open(video)
% writeVideo(video,F)
% close (video)

figure
subplot(211)
stairs(t,u_cl(:,1),'k','linewidth',1.5); axis([0 t(end) -0.2 0.8])
ylabel('v (rad/s)')
grid on
subplot(212)
stairs(t,u_cl(:,2),'r','linewidth',1.5); %axis([0 t(end) -0.85 0.85])
xlabel('time (seconds)')
ylabel('\omega (rad/s)')
grid on