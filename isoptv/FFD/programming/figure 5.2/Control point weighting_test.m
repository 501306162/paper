clear all
close all

% the Control point weighting scheme ,
% seen in paper<Non-rigid Registration Using Free-form Deformations> by Loren Arthur Schwarz
% with the test of the red point in figure 5.2 
[e,f,g,h]=dB(1);
[a,b,c,d]=BB(0);
w=[e,f,g,h];
m=[a,b,c,d];
for i=1:length(w)
    for j=1:length(m)
    q(j,i)=w(i)*m(j);
    end
end
q