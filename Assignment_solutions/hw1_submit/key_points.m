
function exterma=key_points(I)
%I=imread('Image');
I=double(rgb2gray(I));
%I = I(1:end-70,1:end-70);
I = halveSize(I);
I=I/max(max(I));  % image should be in [0 1]
[M,N] = size(I) ;

S=3 ;
omin=-1 ; % first octave -1 mmeans I should be doublsized for first octave
O=floor(log2(min(M,N)))-omin-4 ; % Up to 16x16 images

sigma0=1.6*2^(1/S) ; 
sigman=0.5 
thresh = 0.006;
r = 10 ;

GS = gaussianss(I,O,S,omin,-1,S+1,sigma0) ;
%*************************************************************
%calculate DOG images
for o=1:GS.O %all the octaves
  [M,N,SS] = size(GS.octave{o}) ;
  DOG.octave{o} = zeros(M,N,SS-1) ;
  for s=1:SS-1
    DOG.octave{o}(:,:,s) = ...
        GS.octave{o}(:,:,s+1) - GS.octave{o}(:,:,s) ;
  end
end

%finding key points
exterma=zeros(2,1);
for o=1:GS.O
    for s=2:SS-2
        sig=1.6*2^(o-1)*(2^(1/S))^s;
        current_DOG=DOG.octave{o}(:,:,s);
        down_DOG=DOG.octave{o}(:,:,s-1);
        up_DOG=DOG.octave{o}(:,:,s+1);
 	    extr = search_exterm(up_DOG,down_DOG,current_DOG  ) ;%find exremum
       if extr(1,1)
        extr=localize_eliminate(extr,up_DOG,down_DOG,current_DOG ,thresh,r);
        if extr(1,1)
        extr=2^(o-1+GS.omin) *extr; %stor key points
        exterma=[exterma extr];
        end
       end
    end
end

imshow(I,[])
hold on
plot(exterma(2,:),exterma(1,:),'r+','LineWidth',2)
hold on
exterma = fix(exterma);
[~,Gdir] = imgradient(I);
for i = 2:length(exterma)
    point = [exterma(1,i),exterma(2,i)];
    Pt{i}= estimate_scale_and_orientation(I,Gdir,point,1);
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


end

%**************************************************************************
function SS = gaussianss(I,O,S,omin,smin,smax,sigma0)
%smax--> maximum scale(here it is 4)
%smin--> minimum scale(here is -1...image shold  be double sized)
%omin--> first octave(here is -1)
%1.6 as sigma0 is considered for omin=-1


% Scale multiplicative step
k = 2^(1/S) ;


dsigma0 = sigma0 * sqrt(1 - 1/k^2) ; % Scale step factor
sigman  = 0.5 ;                      % Nominal smoothing of the image

% Scale space structure
SS.O          = O ;
SS.S          = S ;
SS.sigma0     = sigma0 ;
SS.omin       = omin ;
SS.smin       = smin ;
SS.smax       = smax ;

if omin < 0
	for o=1:-omin
		I = doubleSize(I) ;
	end
elseif omin > 0
	for o=1:omin
		I = halveSize(I) ;
	end
end

[M,N] = size(I) ;

% Index offset
so = -smin+1 ;


if(sigma0 * 2^omin * k^smin < sigman)
	warning('The nominal smoothing exceeds the lowest level of the scale space.') ;
end

SS.octave{1} = zeros(M,N,smax-smin+1) ;% we have 6 scale in each octave
SS.octave{1}(:,:,1)  = gauss_filter(I,sqrt((sigma0*k^smin)^2 ...
    - (sigman/2^omin)^2));

for s=smin+1:smax

	dsigma = k^s * dsigma0 ;% smooth Image in prevous scale and just use dsigma
	SS.octave{1}(:,:,s +so) =gauss_filter...
        (squeeze(SS.octave{1}(:,:,s-1 +so)), dsigma);
	%HOSSEIN	imsmooth(squeeze(SS.octave{1}(:,:,s-1 +so)), dsigma )  ;
end


% --------------------------------------------------------------------
%                                                        Other octaves
% --------------------------------------------------------------------

for o=2:O  
	sbest = min(smin + S, smax) ;
	TMP = halveSize(squeeze(SS.octave{o-1}(:,:,sbest+so))) ;
	target_sigma = sigma0 * k^smin ;
	  prev_sigma = sigma0 * k^(sbest - S) ;
            
	if (target_sigma > prev_sigma)
          TMP =gauss_filter(TMP, sqrt(target_sigma^2 - prev_sigma^2));                            
    end
    
    
	[M,N] = size(TMP) ;
	
	SS.octave{o} = zeros(M,N,smax-smin+1) ;
	SS.octave{o}(:,:,1) = TMP ;

	for s=smin+1:smax
		% The other levels are determined as above for the first octave.		
		dsigma = k^s * dsigma0 ;
		SS.octave{o}(:,:,s +so) =gauss_filter(squeeze(SS.octave{o}...
            (:,:,s-1 +so)), dsigma);
	end
	
end
end

% -------------------------------------------------------------------------
%                                                       Auxiliary functions
% -------------------------------------------------------------------------
function J = doubleSize(I)
[M,N]=size(I) ;
J = zeros(2*M,2*N) ;
J(1:2:end,1:2:end) = I ;
J(2:2:end-1,2:2:end-1) = ...
	0.25*I(1:end-1,1:end-1) + ...
	0.25*I(2:end,1:end-1) + ...
	0.25*I(1:end-1,2:end) + ...
	0.25*I(2:end,2:end) ;
J(2:2:end-1,1:2:end) = ...
	0.5*I(1:end-1,:) + ...
    0.5*I(2:end,:) ;
J(1:2:end,2:2:end-1) = ...
	0.5*I(:,1:end-1) + ...
    0.5*I(:,2:end) ;
end

function J = halveSize(I)
J=I(1:2:end,1:2:end) ;
end
%*************************************************************************

%**************************************************************************
function im=gauss_filter(image,sigma)
G = fspecial('gaussian',[5 5],sigma);
im=imfilter(image,G,'same');
end
%**************************************************************************

%**************************************************************************
function [points2]=localize_eliminate(points,up,down,curr,thr,r)
points2=zeros(2,1);
t=1;
for i=1:size(points,2)
    x=points(1,i);
    y=points(2,i);
    fxx= curr(x-1,y)+curr(x+1,y)-2*curr(x,y);   % double derivate in x direction
    fyy= curr(x,y-1)+curr(x,y+1)-2*curr(x,y);   % double derivate in y direction
    fsigmasigma=up(x,y)+down(x,y)-2*curr(x,y);   % double derivate in sigma direction
    
    
    fxsigma=((up(x+1,y)-down(x+1,y))-(up(x-1,y)-down(x-1,y)))/4;%derivate in x and sigma direction
    fysigma=((up(x,y+1)-down(x,y+1))-(up(x,y-1)-down(x,y-1)))/4;%derivate in y and sigma direction
    fxy= curr(x-1,y-1)+curr(x+1,y+1)-curr(x-1,y+1)-curr(x+1,y-1); %derivate inx and y direction
    
    fx=curr(x,y)-curr(x-1,y);%derivate in x direction
    fy=curr(x,y)-curr(x,y-1);%derivate in y direction
    fsigma=(up(x,y)-down(x,y))/2;%derivate in sigma direction
    
    %localization using Teilor seri
    A=[fsigmasigma fxsigma fysigma;fxsigma fxx fxy;fysigma fxy fyy];
    X=-inv(A)*([fsigma fx fy]');
    
    
    x_hat=X(2);
    y_hat=X(3);
    if abs(x_hat)<4 && abs(y_hat)<4 %ignor the ofsets > 4
       px=round(x+x_hat);
       py=round(y+y_hat);
    else
        px=x;
        py=y;
        %[px py]
    end
    
    
    
    D_hat=curr(px,py)+([fsigma fx fy]*X)/2;
    if abs(D_hat)>thr%% filter some low contrast points
        if (fxx+fyy)^2/(fxx*fyy-fxy^2)<(r+1)^2/r % remove edge points
            points2(1,t)=px;points2(2,t)=py;
            t=t+1;
        end
    end
    
    
end
               
end
%**************************************************************************

%**************************************************************************
function indx=search_exterm(up,down,im)
[m n]=size(im);
t=1;
thr=.004;
indx=[0;0];

    for i=2:m-1
        for j=2:n-1
            
        if im(i,j)> thr 
        window(1:3,1:3)=down(i-1:i+1,j-1:j+1);
        window(4:6,1:3)=im(i-1:i+1,j-1:j+1);
        window(7:9,1:3)=up(i-1:i+1,j-1:j+1);
        window(5,2)=-100;
        if im(i,j)>max(max(window))
            indx(:,t)=[i j]';
            t=t+1;
        end
        end
        
        if  im(i,j)<-thr
        window(1:3,1:3)=down(i-1:i+1,j-1:j+1);
        window(4:6,1:3)=im(i-1:i+1,j-1:j+1);
        window(7:9,1:3)=up(i-1:i+1,j-1:j+1);
        window(5,2)=100;
        if im(i,j)<min(min(window))
            indx(:,t)=[i j]';
            t=t+1;
        end
        end
        end
        
    end

end
%**************************************************************************
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
