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
                
                target_positions = hg_target_points([1,2], :) ./ hg_target_points(3, :);
                
                %                 figure(100); imshow(rectified_img); pause(0.001);
                
                for x = 1 : rs_x_size
                    pos = target_positions(:, x);
                    for ch = 1 : rs_channels
                        if (pos(2) <= rs_y_size && pos(2) >= 1 && pos(1) <= rs_x_size && pos(1) >= 1)          
                            
                            % dest_intensity_val = rs_img( round(pos(2)), round(pos(1)), ch);   % this method creates clutter
                            
                            % the following implements bilinear kernel
                            
                            lb = floor(pos);
                            ub = ceil(pos);
                            
                            c1 = [lb(1),  lb(2)];
                            c2 = [lb(1), ub(2)];
                            c3 = [ub(1),  lb(2)];
                            c4 = [ub(1), ub(2)];
                            
                            g1 = ImageRectification.BilinearKernel (pos,  c1);
                            g2 = ImageRectification.BilinearKernel (pos,  c2);
                            g3 = ImageRectification.BilinearKernel (pos,  c3);
                            g4 = ImageRectification.BilinearKernel (pos,  c4);
                            
                            dest_intensity_val = g1 * rs_img(c1(2), c1(1), ch) + g2 * rs_img(c2(2), c2(1), ch) + g3 * rs_img(c3(2), c3(1), ch) + g4 * rs_img(c4(2), c4(1), ch);
                        else
                            dest_intensity_val = rs_img(1, 1, ch) * 0;
                        end
                        rectified_img(y, x, ch) = dest_intensity_val;
                    end
                end
                
            end
            
%             figure(100000); imshow(rectified_img); pause(0.001);
            
        end
        
    end
    
    
    
    methods (Access = private, Static = true)
                
        
        function g = BilinearKernel (uv1, uv2)
            g1 = max(0,   1 - abs(uv1(1) - uv2(1)));
            g2 = max(0,   1 - abs(uv1(2) - uv2(2)));
            g = g1 * g2;
        end
    end
    
end

