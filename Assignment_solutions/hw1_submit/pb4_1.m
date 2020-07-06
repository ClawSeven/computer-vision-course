clear;
%Harris corner detector

image = imread("pb3_1.jpg");
img_ = double(rgb2gray(image));
%img = img_;

%rotation 
%img = rot90(img,1);

%translation 其实就是采样
%img = img_(1:end-200,1:end-200);

%scale of the image
img = halveSize(img_);

[row,col]=size(img);

%gaussian filter parameter
sigma = 1.6;
%N = 3; %kernel size is (2N+1)*(2N+1)
N_row = fix(6*sigma);
gausFilter = fspecial('gaussian',[N_row N_row],sigma);

%%applying sobel edge detector in the horizontal direction
fx = [-1 0 1;-2 0 2;-1 0 1];
Ix = filter2(fx,img);
% applying sobel edge detector in the vertical direction
fy = fx';
Iy = filter2(fy,img); 

Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

% optionaly, gaussion filter
Ix2 = filter2(gausFilter,Ix2);
Iy2 = filter2(gausFilter,Iy2);
Ixy = filter2(gausFilter,Ixy);
R = zeros(row,col);

for i = 1:row
    for j = 1:col
        M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)]; 
        R(i,j) = det(M)-0.05*(trace(M))^2;
    end
end

allmax = max(max(R));
R = R/allmax;
threshold = 0.1;
R(R<threshold)=0;

points = imregionalmax(R);
%points = nlfilter(R, [3 3], @(x) all(x(5)> x([1:4 6:9])) );
[posr, posc] = find(points == 1);
imshow(img,[]);
hold on;
plot(posc,posr,'r.');
hold on;

[~,Gdir] = imgradient(img);



% ptx = posr(10);
% pty = posc(10);
% point = [ptx,pty];
%[best_scale,estimate_response,main_orient]= estimate_scale_and_orientation(img,Gdir,point,1);

              %corner point struct
for i = 1:length(posr)
    point = [posr(i),posc(i)];
    Pt{i}= estimate_scale_and_orientation(img,Gdir,point,1);
    x = Pt{i}.x-Pt{i}.scale;
    y = Pt{i}.y-Pt{i}.scale;
    w = 2*Pt{i}.scale;
    h = 2*Pt{i}.scale;
    
    orientx = Pt{i}.x;
    orienty = Pt{i}.y;
    u = 100*cos(Pt{i}.orient/180*pi);
    v = 100*sin(Pt{i}.orient/180*pi);
    pos = [y x w h];   %%因为图片是旋转了的，所以x位置和y位置要转一下
    
    rectangle('Position',pos,'Curvature',[1 1]);
    hold on;
    
    quiver(orienty,orientx,v,u);
    hold on;
end




function Pt = estimate_scale_and_orientation(img,Gdir,point,dpixel)

    [X,Y] = size(img);
    ptx = point(1);
    pty = point(2);
    max_radius = min(min(X-ptx,ptx-1),min(Y-pty,pty-1));
    scale_range = 1:dpixel:max_radius;
    sigma_range = 0.3*(scale_range-1)+0.8;
    estimate_response = zeros(1,length(scale_range));
    
    for i = 1:length(scale_range)
        patch = img(ptx-scale_range(i):ptx+scale_range(i),pty-scale_range(i):pty+scale_range(i));
        LoGfilter = fspecial('log',[2*scale_range(i)+1,2*scale_range(i)+1],sigma_range(i));
        estimate_response(i) = sum(sum(LoGfilter.*patch));
    end
    
    %argmax, obtain the best scale
    [ ~, idx ] = max(estimate_response);
    best_scale = scale_range(idx);
    
    %calculate the orientation
    orient_patch = Gdir(ptx-best_scale:ptx+best_scale,pty-best_scale:pty+best_scale);
    all_orient = Obtain_orient(orient_patch);  %locate the orientation to certain orientation part
    %statistic of the orientation
    
    hist = tabulate(all_orient(:));
    [~,hist_idx]=max(hist(:,2));
    main_orient = hist(hist_idx,1);
    
    Pt.x = ptx;
    Pt.y = pty;
    Pt.scale = best_scale;
    Pt.orient = main_orient;
    
end

function dirLoc = Obtain_orient(dir)
    if dir>0
        r = fix((dir+5)/10);
        dirLoc = r*10;
    else
        r = fix((dir-5)/10);
        dirLoc = r*10;
    end
end

function J = halveSize(I)
J=I(1:2:end,1:2:end) ;
end









