% function square()
% figure,hold on
% plot(-1:1,-1:1,'-k');
% set(gca,'color',[0.1,0.1,0.1]);
% rectangle('position',[-0.5,-0.5,1,1],'edgecolor','w','facecolor','w')
% axis square 
% hold off
% end


m=256;
n=256;
I=zeros(m,n);

imshow(I,'border','tight','initialmagnification','fit')
axis normal
set (gcf,'Position',[1,1,256,256]);

hold on
imshow(I)

% I=I+rectangle('position',[256/4,256/4,128,128],'edgecolor','w','facecolor','w') ;   % Ìî³äÉ«
fill([256/4,256*3/4,256*3/4,256/4],[256/4,256/4,256*3/4,256*3/4],'w') 
 hold off
 
saveas(gca,'square.bmp','bmp');