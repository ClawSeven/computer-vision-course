function y = Homo_tran(x,v)
% a series of points
q = v * [x; ones(1, size(x,2))];
p = q(3,:);
y = [q(1,:)./p; q(2,:)./p];
end