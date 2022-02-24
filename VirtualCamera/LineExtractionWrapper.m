classdef LineExtractionWrapper<handle
    properties (Access = protected)
        num_windows = 0;
        rejection_identifier = 0;
        inlier_threshold = 0;
        doPlotResults = 0;
        addedNoiseLevels = 0;
        rs_image = 0;
        template = 0;
        K = 0;
    end
    properties (Constant)
        point_size = 32;
        font_size = 22;
        line_width = 2;
    end
    methods (Access = public)
        function [this] = LineExtractionWrapper(num_windows, rejection_identifier, ...
                                                                inlier_threshold, addedNoiseLevels, doPlotResults,... 
                                                                                                                        rs_image, template, K)
            this.num_windows = num_windows;
            this.rejection_identifier = rejection_identifier;
            this.inlier_threshold = inlier_threshold;
            this.doPlotResults = doPlotResults;
            this.addedNoiseLevels = addedNoiseLevels;
%             this.rs_image = rs_image;
%             this.template = template;
%             this.K = K;
        end
        
        function [refined_lines] = extract_refined_lines(this, retinal_matches, template_matches)
            [horizon] = HorizontalLineExtraction(this.num_windows, this.rejection_identifier, this.inlier_threshold);
            [scanline_ids] = horizon.extract_line(retinal_matches, ones(size(retinal_matches,1),1));
            D1 = retinal_matches(scanline_ids~=this.rejection_identifier,1:2);
            D2 = template_matches(scanline_ids~=this.rejection_identifier,1:2);
            ids = scanline_ids(scanline_ids~=this.rejection_identifier,1);
            D1n = D1+(randn(size(D1,1),2)*this.addedNoiseLevels);
            if(this.doPlotResults == 1)
                figure;
                set(gca,'FontSize',this.font_size);
                subplot(1,2,1);
                scatter(D1(:,1),D1(:,2),this.point_size, 'b+'); hold on;
                scatter(D1n(:,1),D1n(:,2),this.point_size, 'r*'); hold on;
                unique_scaline = unique(ids);
                min_retinal = min(D1(:,1));
                max_retinal = max(D1(:,1));
                for ii = 1:size(unique_scaline,1)
                    plot([min_retinal max_retinal],[unique_scaline(ii,1) unique_scaline(ii,1)],'LineWidth',this.line_width);
                end
                legend('Without noise','With noise');
                xlabel('X-axis');
                ylabel('Y-axis');
                grid on;
                hold off; drawnow;
                title('Retinal keypoints + scanlines');
                subplot(1,2,2);
                scatter(D2(:,1),D2(:,2), this.point_size, 'r+'); drawnow;
                title('Template keypoints');
                xlabel('X-axis');
                ylabel('Y-axis');
                set(findall(gcf,'-property','FontSize'),'FontSize',this.font_size);
            end
            
            refined_lines.retinal_matches = D1n;
            refined_lines.template_matches = D2;
            refined_lines.ids = ids;%this.unnormalise_ids(ids);
        end
        
        function [id_ref] = unnormalise_ids(this, ids)
            id_ref = zeros(size(ids,1),1);
            for ii = 1:size(ids,1)
                p = [0;ids(ii);1];
                pU = this.K*p;
                id_ref(ii) = round(pU(2));
            end
        end
        
        function draw_matches(this, refined_lines)
            figure; ax = axes;
            D1 = refined_lines.retinal_matches;
            D2 = refined_lines.template_matches;
            D1p = [D1 ones(size(D1,1),1)];
            D1un = this.K*D1p.';
            D1un(3,:) = [];
            showMatchedFeatures(this.rs_image,this.template,D1un.', D2 ,'montage','Parent',ax);
        end
    end
end