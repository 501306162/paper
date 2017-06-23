% % function circle()
m=256;
n=256;
I=zeros(m,n);

imshow(I,'border','tight','initialmagnification','fit')
axis normal
set (gcf,'Position',[1,1,256,256]);

hold on
imshow(I)
alpha=0:pi/50:2*pi;%角度[0,2*pi] 
R=256/4;%半径 
x=R*cos(alpha); 
y=R*sin(alpha); 
    
fill(x+128,y+128,'w')        % 填充色
 hold off
saveas(gca,'circle.bmp','bmp');

% end
% 
% clear;clc
% m=256;
% n=256;
% I=zeros(m,n);
% 
% hold on ;
% imshow(I);
% alpha=0:pi/50:2*pi;%角度[0,2*pi] 
% R=0.5;%半径 
% x=R*cos(alpha); 
% y=R*sin(alpha); 
% 
% fill(x,y,'w')
% axis square
% hold off;
% set(gcf,'color','k');
% % set(gca,'units','pixels','Visible','off');
% q=get(gca,'position');
% q(1)=0;%设置左边距离值为零
% q(2)=0;%设置右边距离值为零
% set(gca,'position',q);
% frame=getframe(gcf,[1,1,n,m]);%
% im=frame2im(frame);
% Ibw=im2bw(im,graythresh(im));
% figure,imshow(Ibw);
% imwrite(Ibw,'a.jpg','jpg');%可以修改保存的格式