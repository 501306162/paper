
% clc,clear all,close all;

%% Read two images  
Im=im2double(imread('images/statue-rio-deformed.png'));
Is=im2double(imread('images/statue-rio.png'));

ref = Is;
frame = Im ;

[B1,A1] = demosimp(frame,ref);

upref1 = movepixels(frame,B1,A1);

figure;
subplot(1,3,1),imshow(ref,[]);title('static');
subplot(1,3,2),imshow(frame,[]);title('moving');
subplot(1,3,3),imshow(upref1,[]);title('registered');

I1 = ref - frame;
I2 = ref - upref1;

figure;
subplot(1,2,1),imshow((I1),[]);title('‘≠≤Ó');
subplot(1,2,2),imshow((I2),[]);title('œ÷≤Ó');

mse=myMSE( ref,upref1 )
tss=now();
names=strcat(num2str(hour(tss)),'_',num2str(minute(tss)));
namex=strcat('Yp_',names,'.mat');
fprintf(strcat(namex,'\n'))
save(namex,'mse','upref1');

