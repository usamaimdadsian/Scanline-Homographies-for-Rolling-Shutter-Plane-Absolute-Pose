classdef YizhenRS < handle
    
    properties(Access=private)
        
        isRSImageAvailable = false;
        isTemplateImageAvailable = false;
        isKeyPointsSet = false;
        isIntrinsicsSet = false;
        
        keypoints_rs = []; % expecting RS keypoints
        keypoints_template = []; % expecting template keypoints        
        
        image_rs = []; % to load the RS image
        image_template = []; % to load the template image
        
        intrinsics = [];
        
        rectifiedImage = [];
        landmarkPoses = [];
        landmark_homography = [];
        additional_landmarks_RS = [];
    end
    
    methods(Access=public)
        
        function [isDone] = rectify(obj)
            assert((obj.isRSImageAvailable && obj.isTemplateImageAvailable && obj.isKeyPointsSet && obj.isIntrinsicsSet )...
                ,'Either the keypoints or the RS/template images or intrinsics have not been set. To check usage, hit: "help YizhenRS" ');
            
            [obj.rectifiedImage, obj.landmarkPoses, obj.landmark_homography] = TwoViewRectification(obj.keypoints_rs, obj.keypoints_template,...
                                                                            obj.image_rs, obj.image_template, obj.intrinsics, obj.additional_landmarks_RS);            
                                                                        
            isDone=true;
        end
        
        
        function [obj] = YizhenRS(keypoints_template_, keypoints_rs_)            
            is_kps1_set = false; is_kps2_set = false;
            
            if exist('keypoints_rs_','var')
                obj.keypoints_rs = keypoints_rs_;
                is_kps1_set = true;
            end
            
            if exist('keypoints_template_','var')
                obj.keypoints_template = keypoints_template_;
                is_kps2_set = true;
            end
            
            if(is_kps1_set && is_kps2_set)
                obj.isKeyPointsSet = true;
            end            
        end

        function [obj] = set_template_image(obj, image)
            obj.image_template = image;
            obj.isTemplateImageAvailable = true;
        end
        
        function [obj] = set_RS_image(obj, image)
            obj.image_rs = image;
            obj.isRSImageAvailable = true;
        end
        
        function [obj] = set_intrinsics(obj, K)
            obj.intrinsics = K;
            obj.isIntrinsicsSet = true;
        end
        
        function [obj] = set_additional_landmarks(obj, additional_landmarks_RS_)
            obj.additional_landmarks_RS = additional_landmarks_RS_;
        end
        
        function [image] = get_rectified_image(obj)
            image = obj.rectifiedImage;
        end
        
        function [landmarks] = get_landmark_poses(obj)
            landmarks = obj.landmarkPoses;
        end
        
        function [landmarkHomography] = get_landmark_homographies(obj)
            landmarkHomography = obj.landmark_homography;
        end        

        
    end
end