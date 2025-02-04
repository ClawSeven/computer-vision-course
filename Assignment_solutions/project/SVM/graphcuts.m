function [Ncut] = graphcuts(I,pad,MAXVAL)
% function [Ncut] = graphcuts(I)
% Input: I image
% pad: spatial connectivity; eg. 3
% MAXVAL: maximum image value
% Output: Ncut: Binary map 0 or 1 corresponding to image segmentation
I = double(I); [H,W] = size(I);
% Find weights between nodes I1 and I2, w = exp(a*abs(I1-I2));
% Set a to have a weight of 0.01 for diff = MAXVAL
a = log(0.01)/MAXVAL; x = [0:MAXVAL/100:MAXVAL]'; y = exp(a*x);
figure;plot(x,y);xlabel('intensity diff');ylabel('weights'); title('weights')
ws = 2*pad + 1;
if(ws <= 3)
    ws = 3;
end
%Build the weight matrix
disp('Building Weight Matrix'); close all; tic
WM = zeros(H*W,H*W); countWM = 0;
for kk = 1:W
for jj = 1:H
mask = logical(zeros(H,W));
cs = kk-pad; ce = kk+pad; rs = jj-pad; re = jj+pad;
if(cs<1) 
    cs = 1; 
end;
if(ce>W)
    ce = W; 
end;
if(rs<1) 
    rs = 1;
end;
if(re>H) 
    re = H; 
end;
mask(rs:re,cs:ce) = 1;
idx = find(mask==1);
p = abs(I(idx) - I(jj,kk)); p = exp(a*p);
countWM = countWM + 1; WM(countWM,idx) = p(:)';
end
end
ttime = toc; disp(sprintf('Time for generating weight matrix = %f',ttime)); clear countWM
% Weight between a node and iteself is 0
for jj = 1:H*W 
    WM(jj,jj) = 0; 
end; 
WM = sparse(WM);
% Shi and Malik Algorithm: second smallest eigen vector
disp('Finding Eigen Vector');
d = sum(WM,2); D = diag(d); tic
B = (D-WM); B = (B+B')/2; OPTS.disp = 0;
[v,d,flag] = eigs(B,D,2,'SA',OPTS); ttime = toc;
disp(sprintf('Time for finding eigen vector = %f',ttime)); clear OPTS
y = v(:,2);
Ncut = reshape(y,H,W);
Ncut = Ncut > 0;