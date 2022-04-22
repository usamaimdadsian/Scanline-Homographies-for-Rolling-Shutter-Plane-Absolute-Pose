classdef ImageRectification
   
    methods (Static = true)
        % Method to rectify images using supplied planar homographies per row of image
        % PARAMETERS:
        %  rs_img - rolling shutter image
        %  plane_homographies - cell(s) containing the Homography matrices
        %  anchor_idx - the anchor (row of image) against which the
        %                           rectification is needed
        %
        % RETURNS:
        %  rectified_img - rectifiec image of the same class as the input
        %  points2d - [N x 2] array of points containing just the rectified
        %                       pixel positions (helps in debugging)
        function [rectified_img, points2d] = rectifyRSImage (rs_img, plane_homographies, anchor_idx)

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
            H0 = plane_homographies{anchor_idx};
            fprintf(1, 'The anchor scanline index is %d\n', anchor_idx);

            x_coords = 1 : rs_x_size;
            all_ones = ones(1, rs_x_size);


            % this accumlates interpolation weights.
            weight_matrix = zeros (rs_y_size, rs_x_size);

            for y = 1 : rs_y_size

                rectH = H0 * inv(plane_homographies{y}) ;
                hg_scanline_points = [x_coords;  y*all_ones;  all_ones];
                hg_target_points = rectH * hg_scanline_points;
                target_positions = hg_target_points([1,2], :) ./ hg_target_points(3, :);
                                
                for x = 1 : rs_x_size
                    pos = target_positions(:, x);
                    if (pos(2) <= rs_y_size && pos(2) >= 1 && pos(1) <= rs_x_size && pos(1) >= 1)
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

                        weight_matrix(c1(2), c1(1)) = weight_matrix(c1(2), c1(1)) + g1;
                        weight_matrix(c2(2), c2(1)) = weight_matrix(c2(2), c2(1)) + g2;
                        weight_matrix(c3(2), c3(1)) = weight_matrix(c3(2), c3(1)) + g3;
                        weight_matrix(c4(2), c4(1)) = weight_matrix(c4(2), c4(1)) + g4;
                    end
                end

            end


            % initialize the rectified image size
            rectified_img = rs_img * 0;

            points2d = zeros(size(rs_img,1)*size(rs_img,2),2);
            start = 1; endI = rs_x_size;

            for y = 1 : rs_y_size
                rectH = H0 * inv(plane_homographies{y}) ;

                hg_scanline_points = [x_coords;  y*all_ones;  all_ones];

                hg_target_points = rectH * hg_scanline_points;

                target_positions = hg_target_points([1,2], :) ./ hg_target_points(3, :);
                
                for x = 1 : rs_x_size
                    pos = target_positions(:, x);
                    for ch = 1 : rs_channels
                        if (pos(2) <= rs_y_size && pos(2) >= 1 && pos(1) <= rs_x_size && pos(1) >= 1)
                            dest_intensity_val = rs_img(y, x, ch);
                            % we need to set pixel intensity to pos(2), pos(1)
                            % which is fractional                            
                            % rectified_img(pos(2), pos(1), ch) = dest_intensity_val;

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

                            w1 = weight_matrix(c1(2), c1(1));
                            w2 = weight_matrix(c2(2), c2(1));
                            w3 = weight_matrix(c3(2), c3(1));                            
                            w4 = weight_matrix(c4(2), c4(1));

                            if (true)
                                if (w1 > 0)
                                    g1 = g1/w1;
                                end
                                if (w2 > 0)
                                    g2 = g2/w2;
                                end
                                if (w3 > 0)
                                    g3 = g3/w3;
                                end
                                if (w4 > 0)
                                    g4 = g4/w4;
                                end
                            end
                            
                            rectified_img(c1(2), c1(1), ch) = rectified_img(c1(2), c1(1), ch) + g1 * dest_intensity_val;
                            rectified_img(c2(2), c2(1), ch) = rectified_img(c2(2), c2(1), ch) + g2 * dest_intensity_val;
                            rectified_img(c3(2), c3(1), ch) = rectified_img(c3(2), c3(1), ch) + g3 * dest_intensity_val;                            
                            rectified_img(c4(2), c4(1), ch) = rectified_img(c4(2), c4(1), ch) + g4 * dest_intensity_val;
                            
                        end
                    end
                end

            end

        end



        function rectified_img = ImageInterpolationByNeighborhood (rectified_img)

            [y_size, x_size, channels] = size(rectified_img);

            % use interpolation to fill holes
            for y = 2 : y_size-1
                for x = 2 : x_size-1
                    for ch = 1 : channels
                        dest_intensity_val = rectified_img(y, x, ch);
                        if (~dest_intensity_val)

                            % the following implements average of the neighbourhood
                            pos = [x, y];

                            c1 = [x-1,  y];
                            c2 = [x+1, y];
                            c3 = [x,  y-1];
                            c4 = [x, y+1];

                            p1 = rectified_img(c1(2), c1(1), ch);
                            p2 = rectified_img(c2(2), c2(1), ch);
                            p3 = rectified_img(c3(2), c3(1), ch);
                            p4 = rectified_img(c4(2), c4(1), ch);

                            nocc = nnz([p1, p2, p3, p4]);
                            if nocc > 2
                                ss = 1/nocc;
                                dest_intensity_val = ss * p1 + ss * p2 + ss * p3 + ss * p4;
                            end

                            rectified_img(y, x, ch) = dest_intensity_val;

                        end
                    end
                end
            end

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

