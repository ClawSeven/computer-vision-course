function [res,index]= Warp_Image(img1,img2,H_Matrix)
%%%����һ��ͼimg2���Ͱ����任����һ��ͼ�ĵ�Ӧ�Ծ���H_Matrix,�õ��任֮���ͼ
%%%Ȼ���img1��img2�ϳ�һ��
%%%ʹ�õ���reverse warp
%%%index ���ص�����Ϊframe/ reference������ͼԭ����1,1���ڱ任��Ĵ�ͼ�е�λ��index [��;��]����[x,y]

[row2, col2] = size(img2(:,:,1));
[row1, col1] = size(img1(:,:,1));
%%%�߽ǵ�λ��
corner_P = [1,col2,1,col2;1,1,row2,row2];
change_P = Homo_tran(corner_P,H_Matrix);
change_P = round(change_P);

x_min = min([min(change_P(1,:)),1]);
y_min = min([min(change_P(2,:)),1]);
x_max = max([max(change_P(1,:)),col1]);
y_max = max([max(change_P(2,:)),row1]);

[X,Y] = meshgrid(x_min:x_max,y_min:y_max);
[row3,col3] = size(X);

%%%%ʹ��reshapeһ�μ��������
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
P = [X;Y];
P_tran = Homo_tran(P,inv(H_Matrix));
tran_X = reshape(P_tran(1,:),row3,col3);
tran_Y = reshape(P_tran(2,:),row3,col3);
index = [2-x_min;2-y_min];
res = zeros(row3,col3,3);
for color = 1:3
    output = interp2(double(img2(:,:,color)),tran_X,tran_Y);
    output = round(output);
    nan_pos = isnan(output);
    output(nan_pos) = 0;%%%img2�任֮���ͼ��
    output2 = zeros(row3,col3);%%img1�ڸ�ͼ�е�λ�ã�
    output2((2-y_min):(1-y_min+row1),(2-x_min):(1-x_min+col1))=double(img1(:,:,color));
    pos = (output2~=0);
    output(pos)=0;
    %pos = (output~=0);
    %output2(pos)=0;
    res(:,:,color) = imadd(output,output2);
end
res = uint8(res);
end