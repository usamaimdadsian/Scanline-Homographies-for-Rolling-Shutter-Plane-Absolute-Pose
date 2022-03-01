classdef MotionSimulator < handle

    properties (Access = public)
        
        mPoseCell = {};
        
    end
    
    
    
    methods (Access = public)
        
        function obj = MotionSimulator (InitialPose)
            obj.init_R = InitialPose(1:3, 1:3);
            obj.init_t = InitialPose(1:3, 4);
            obj.mPoseCell{1} = [obj.init_R, obj.init_t];
        end
        
    end
    
    
    methods (Access = public)
                
        function simPoseCell = SimuStaticMotionGaussianNoise (this,  sigma_rot, sigma_pos, numPoses)
        % static at the original pose with Gaussian noise
            R = this.init_R;
            t = this.init_t;
            n_exist_num_poses = length(this.mPoseCell);
            for ii = (n_exist_num_poses+1):(n_exist_num_poses+numPoses)            
                R = R * internal_packages.SO3.Exp(sigma_rot *randn(3, 1));
                t = t + sigma_pos *randn(3, 1);
                this.mPoseCell{ii} = [R, t];
            end
            this.init_R = R;   this.init_t = t;
            simPoseCell = this.mPoseCell;
        end
        

        function simPoseCell = SimuStaticMotionRandomWalk (this,  sigma_rot, sigma_pos, numPoses)
        % static at the original pose with random walk
            R = this.init_R;
            t = this.init_t;
            n_exist_num_poses = length(this.mPoseCell);
            for ii = (n_exist_num_poses+1):(n_exist_num_poses+numPoses)
                R = R * internal_packages.SO3.Exp(sigma_rot * (rand(3, 1) - 0.5));
                t = t + sigma_pos * (rand(3, 1) - 0.5);
                this.mPoseCell{ii} = [R, t];
            end
            this.init_R = R;   this.init_t = t;
            simPoseCell = this.mPoseCell;
        end        
        
        
        function simPoseCell = SimuConstVelocityMotionGaussianNoise (this, velocity_rot, velocity_pos, sigma_rot, sigma_pos, numPoses)
        % constant velocity with Gaussian noise
            R = this.init_R;
            t = this.init_t;
            n_exist_num_poses = length(this.mPoseCell);
            for ii = (n_exist_num_poses+1):(n_exist_num_poses+numPoses)
                R = R * internal_packages.SO3.Exp(velocity_rot) * internal_packages.SO3.Exp(sigma_rot * randn(3, 1));
                t = t + velocity_pos + sigma_pos * randn(3, 1);
                this.mPoseCell{ii} = [R, t];
            end
            this.init_R = R;   this.init_t = t;
            simPoseCell = this.mPoseCell;
        end
        

        function simPoseCell = SimuConstVelocityMotionRandomWalk (this, velocity_rot, velocity_pos, sigma_rot, sigma_pos, numPoses)
        % constant velocity with random walk
            R = this.init_R;
            t = this.init_t;
            n_exist_num_poses = length(this.mPoseCell);
            for ii = (n_exist_num_poses+1):(n_exist_num_poses+numPoses)
                R = R * internal_packages.SO3.Exp(velocity_rot + sigma_rot * (rand(3, 1) - 0.5));
                t = t + velocity_pos + sigma_pos * (rand(3, 1) - 0.5);
                this.mPoseCell{ii} = [R, t];
            end
            this.init_R = R;   this.init_t = t;
            simPoseCell = this.mPoseCell;
        end        
        
        
    end

    
    % the very fist pose
    properties (Access = private)
        
        init_R = []
        
        init_t = []
        
    end
        
    
end

