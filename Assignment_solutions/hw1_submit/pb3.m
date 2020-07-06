clear;
image1 = "pb3_1.jpg";
image2 = "pb3_2.jpg";

img1_ = imread(image1);
img2_ = imread(image2);

img1 = double(rgb2gray(img1_));
img2 = double(rgb2gray(img2_));
%[row,col] = size(graycat);

sigma = 10;
N = fix(3*sigma);
N_row = 2*N+1;
gausFilter = fspecial('gaussian',[N_row N_row],sigma);

blur1 = filter2(gausFilter,img1);
blur2 = filter2(gausFilter,img2);

out = blur1+(img2-blur2);
imshow(out,[]);

% blur1 = imfilter(img1,gausFilter,'conv');
% close = graycat-blur;
% imshow(close);
% 
% subplot(1,2,1);imshow(blur);title('Ô¶¿´');
% subplot(1,2,2);imshow(close);title('½ü¿´');




%imshow(cat);
%imshow(dog);

