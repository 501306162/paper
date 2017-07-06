function neighbors_ind=neighborhood_trans(ind,image,K)
% return the neighbors'index of 'ind ' pixel 
N=numel(size(K));
switch (N)
    case 2
       neighbors_ind= neighborhood2d(ind,image,K);
    case 3
       neighbors_ind= neighborhood3d(ind,image,K);   
    otherwise
       fprintf('%d is not include in the neighborhood type!',N); 
end

end
function neighbors_ind=neighborhood2d(ind,image,K)
% return the neighbors'index of 'ind ' pixel
[r,c]=ind2sub(size(image),ind);
subfirst=[K(1),K(2)];
subend=size(image);

% subscript in moving image
neig_sub_begin=max(subfirst,[r,c]-subfirst);
neig_sub_end=min(subend,[r,c]+subfirst);

row_size=neig_sub_end(1)-neig_sub_begin(1)+1;
col_size=neig_sub_end(2)-neig_sub_begin(2)+1;

row=repmat((neig_sub_begin(1):1:neig_sub_end(1)),col_size,1);   % row(:)-->111 222 333
col=repmat((neig_sub_begin(2):1:neig_sub_end(2)),1,row_size);   % col(:)-->123 123 123

% the neighbors include the 'ind' pixel
neighbors_ind=sub2ind(size(image),row(:),col(:));
% the neighbors exclude the 'ind' pixel
neighbors_ind(neighbors_ind==ind)=[];
% % translate the current pixel to the first one
% neighbors_ind=[ind,neighbors_ind(:)];

% % form the cliques
% clique
end

function neighbors_ind=neighborhood3d(ind,image,K)
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