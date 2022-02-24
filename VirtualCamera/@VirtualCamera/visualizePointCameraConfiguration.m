function hFig = visualizePointCameraConfiguration (fid, CameraMatrices, Points3D)

hFig = figure(fid);

if (nargin == 3)
    scatter3(Points3D(1,:), Points3D(2,:), Points3D(3,:), '.');
end

hold on;

for ii = 1 : length(CameraMatrices)
    
    P = CameraMatrices{ii};
    
    if ( abs(det(P(1:3, 1:3))) < 1e-12)
        
        fprintf(2, '\nCamera [%d] lies at infinity. det(P33) = %.18f \n', ii, abs(det(P(1:3, 1:3))));
        fprintf(2, 'Camera at infity is not plotted in the visualization. \n\n');
        
    else
        
        [K, R, C] = VirtualCamera.decomposeCameraMatrix (P);
        plotCamera('Location', C, 'Orientation', R, 'Opacity', 0, 'Color', 'k', 'AxesVisible', true, 'Size', 10);
        
    end
    
end

hold off;
axis equal;
view(3);

pause(0.1);

end