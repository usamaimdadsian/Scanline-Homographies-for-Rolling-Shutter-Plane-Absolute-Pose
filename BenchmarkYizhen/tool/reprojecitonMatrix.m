function  reprojecitonMatrix(K,lines,R,t,w,d,image,color,pathImage)

for i = 1:size(lines,2)
    for j = 1:size(lines{i},2)
        temp = RS_PrejectionPredication([lines{i}(:,j);0],K,R,t,X_(w),d);
        lines{i}(:,j) = temp(1:2,1);
    end
end

IFigure = figure;
%image = rgb2gray(image);
imshow(image,'border','tight','initialmagnification','fit');
n = size(image,1);
m = size(image,2);
set (gcf,'Position',[0,0,m,n]);
axis normal;

figure(IFigure)
hold on
for i = 1:size(lines,2)
    plot(lines{i}(1,:),lines{i}(2,:),'Color',color,'LineStyle','-','LineWidth',3,'Marker','none');
end
hgexport(gcf, pathImage, hgexport('factorystyle'), 'Format', 'jpeg');
close

