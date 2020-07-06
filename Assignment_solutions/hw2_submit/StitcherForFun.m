clear;
image1 = imread("referenced_image.jpg");%%reference image
image2 = imread("transformed_image.jpg");%%need being tranformed image

[mp,fp] = cpselect(image2,image1,'Wait',true);
mp = mp';
fp = fp';
H_Matrix = Homo_solve(mp,fp);

[output,index] = Warp_Image(image1,image2,H_Matrix);
imshow(output,[]);
hold on
plot(fp(1,:)+index(1)-1,fp(2,:)+index(2)-1,'r+','LineWidth',2);
hold on;
img2_tran = Homo_tran(mp,H_Matrix);
plot(img2_tran(1,:)+index(1)-1,img2_tran(2,:)+index(2)-1,'bo','LineWidth',2);

% imshow(image1,[]);
% hold on
% plot(fp(1,:),fp(2,:),'bo','LineWidth',2);
