function ddk=neighborhood_trans(ind,sizes,K)
% return the neighbors'index of 'ind ' pixel 
N=numel(size(K));
switch (N)
    case 2
       ddk= neighborhood2d(ind,sizes,K);
    case 3
       ddk= neighborhood3d(ind,sizes,K);   
    otherwise
       fprintf('%d is not include in the neighborhood type!',N); 
end

end
function ddk=neighborhood2d(ind,sizes,K)
% return the neighbors'index of 'ind ' pixel
[r,c]=ind2sub(sizes,ind);
subfirst=[K(1),K(2)];
subend=sizes;

% subscript in moving image
neig_sub_begin=max(subfirst,[r,c]-subfirst);
neig_sub_end=min(subend,[r,c]+subfirst);

row_size=neig_sub_end(1)-neig_sub_begin(1)+1;
col_size=neig_sub_end(2)-neig_sub_begin(2)+1;

row=repmat((neig_sub_begin(1):1:neig_sub_end(1)),col_size,1);   % row(:)-->111 222 333
col=repmat((neig_sub_begin(2):1:neig_sub_end(2)),1,row_size);   % col(:)-->123 123 123

% the neighbors include the 'ind' pixel
neighbors_ind=sub2ind(sizes,row(:),col(:));
% the neighbors exclude the 'ind' pixel
neighbors_ind(neighbors_ind==ind)=[];

ddk_1=max((1-abs(row-r)./K(1)),0);
ddk_2=max((1-abs(col-c)./K(2)),0);

% ddk=zeros(num_control,1);
ddk=[neighbors_ind,ddk_1.*ddk_2];
% % translate the current pixel to the first one
% neighbors_ind=[ind,neighbors_ind(:)];

% % form the cliques
% clique
end

function ddk=neighborhood3d(ind,image,K)
% return the neighbors'index of 'ind ' pixel

[r,c,s]=ind2sub(size(image),ind);
subfirst=[K(1),K(2),K(3)];
subend=size(image);

% subscript in moving image
neig_sub_begin=max(subfirst,[r,c,s]-subfirst);
neig_sub_end=min(subend,[r,c,s]+subfirst);

row_size=neig_sub_end(1)-neig_sub_begin(1);
col_size=neig_sub_end(2)-neig_sub_begin(2);
slice_size=neig_sub_end(3)-neig_sub_begin(3);

row=repmat((neig_sub_begin(1):1:neig_sub_end(1)),col_size*slice_size,1);            % row(:)-->(111...222...333...)
col=repmat(repmat((neig_sub_begin(2):1:neig_sub_end(2)),row_size,1),1,slice_size);% col(:)-->(111222333...111222333)
slice=repmat((neig_sub_begin(3):1:neig_sub_end(3)),1,row_size*col_size);        % slice(:)-->(123...123...123...)

% the neighbors include the 'ind' pixel
neighbors_ind=sub2ind(size(image),row(:),col(:),slice(:));
% the neighbors exclude the 'ind' pixel
neighbors_ind(neighbors_ind==ind)=[];
end