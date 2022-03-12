function  reprojecitonMatrixCaley(cameraParams,lines, Rot_Rows, trans_Rows,image, color,pathImage)

height = cameraParams.IntrinsicMatrix(3,2) * 2;

for i = 1:size(lines,2)
    for j = 1:size(lines{i},2)
        temp = RSWorld2Image([lines{i}(:,j)' 0], Rot_Rows, trans_Rows, cameraParams);
        lines{i}(:,j) = temp(1,1:2)';
    end
end

IFigure = figure;
image = rgb2gray(image);
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

