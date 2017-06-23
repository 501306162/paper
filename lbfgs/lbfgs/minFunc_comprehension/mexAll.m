% minFunc
% cd ('E:\GitHub\paper\lbfgs\lbfgs\minFunc_comprehension')
% addpath(genpath(pwd))
% mexAll
fprintf('Compiling minFunc files...\n');
cd('minFunc/mex/movepixels');
files=dir('*.c');
for i=1:length(files)
    filename=[files(i).name];
    if(length(filename)<19||(~strcmpi(filename(1:19),'image_interpolation')))
        disp(['compiling : ' filename]);
        mex(filename,'image_interpolation.c');
    end
end
cd('..');
cd('E:/GitHub/paper/lbfgs/lbfgs/minFunc_comprehension');
mex -outdir minFunc/compiled minFunc/mex/mcholC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
mex -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c


