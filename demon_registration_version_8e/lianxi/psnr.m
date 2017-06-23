function S=psnr(I,J)
% PSNR(I,J) returns the peak signal to noise ratio) between I and J %(dB)
% I is the original image and J is a modified version of I. 
% The PSNR value is useful to calculate the distortions on an image.

if (size(I)~=size(J))
   error('Size mismatch!')
end

   [m n] = size(I);
   A=double(I);
   B=double(J);
   sumaDif=0;
   maxI=m*n*max(max(A.^2));
   sumaDif=sum(sum((A-B).^2));   
   if (sumaDif==0) 
      sumaDif=1;
   end
   S=maxI./sumaDif;
   S=10*log10(S);	

end
