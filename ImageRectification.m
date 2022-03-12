classdef ImageRectification
    
    methods (Static = true)
    
        function rectified_img = rectifyRSImage (rs_img, plane_homographies, anchor_idx)

            [rs_y_size, rs_x_size, rs_channels] = size(rs_img);
            
            % verify the rolling-shutter direction. 
            % typically, the rolling shutter direction is from top to bottom
            if (length(plane_homographies) == rs_y_size)
                fprintf(1, 'The rolling shutter direction is from-top-to-bottom\n');
            elseif (length(plane_homographies) == rs_x_size)
                fprintf(2, 'warning@ImageRectification: The rolling shutter direction is from-left-to-right\n');
                fprintf(2, 'This is unusual!\n Press any key to continue if this is really the case\n');
                pause;
            else
                fprintf(2, 'error@ImageRectification: The size of plane_homographies is incorrect!\n'); return;
            end
            
            % choose the anchor scanline
            if ~exist('anchor_idx', 'var')
                anchor_idx = round(rs_y_size/2);
            end            
            invH0 = pinv(plane_homographies{anchor_idx});
            fprintf(1, 'The anchor scanline index is %d\n', anchor_idx);
            
            % initialize the rectified image size
            rectified_img = rs_img * 0;
            
            x_coords = 1 : rs_x_size;
            all_ones = ones(1, rs_x_size);
            
            
            for y = 1 : rs_y_size
                
                rectH = plane_homographies{y} * invH0;
                
                hg_scanline_points = [x_coords;  y*all_ones;  all_ones];
                
                hg_target_points = rectH * hg_scanline_points;
                
                hg_target_points([1,2], :) = hg_target_points([1,2], :) ./ hg_target_points(3, :);
                
                target_positions = round(hg_target_points([1,2], :));
                
%                 figure(100); imshow(rectified_img); pause(0.001);
                
                for x = 1 : rs_x_size
                    pos = target_positions(:, x);
                    for ch = 1 : rs_channels                        
                        if (pos(2) <= rs_y_size && pos(2) >= 1 && pos(1) <= rs_x_size && pos(1) >= 1)
                            dest_intensity_val = rs_img(pos(2), pos(1), ch);
                        else
                            dest_intensity_val = rs_img(1, 1, ch) * 0;
                        end 
                        rectified_img(y, x, ch) = dest_intensity_val;
                    end
                end
                
            end
            
        end
        
    end
    
end

