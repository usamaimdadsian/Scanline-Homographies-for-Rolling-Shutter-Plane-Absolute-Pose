clear all
close all
clc


InitPose = eye(4);



MS = internal_packages.MotionSimulator(InitPose);



sigma_rot = 0.01; sigma_pos = 0.01;
simulation_steps = 9;
poses = MS.SimuStaticMotionGaussianNoise (sigma_rot, sigma_pos, simulation_steps);


sigma_rot = 0.01; sigma_pos = 0.01;
simulation_steps = 10;
poses = MS.SimuStaticMotionRandomWalk (sigma_rot, sigma_pos, simulation_steps);


sigma_rot = 0.01; sigma_pos = 0.01;
simulation_steps = 10;
velocity_rot = [0.1; 0.1; 0.1]; velocity_pos = [0.01; 0.01; 0.01];
poses = MS.SimuConstVelocityMotionGaussianNoise (velocity_rot, velocity_pos, sigma_rot, sigma_pos, simulation_steps);


sigma_rot = 0.01; sigma_pos = 0.01;
simulation_steps = 10;
velocity_rot = [0.1; 0.1; 0.1]; velocity_pos = [0.01; 0.01; 0.01];
poses = MS.SimuConstVelocityMotionRandomWalk (velocity_rot, velocity_pos, sigma_rot, sigma_pos, simulation_steps);


figure(1); hold on;
h1 = plotPoses (poses(1:10))
h2 = plotPoses (poses(10:20))
h3 = plotPoses (poses(20:30))
h4 = plotPoses (poses(30:40))
legend([h1, h2, h3, h4], {'Static Gaussian', 'Static Random walk', 'Velocity Gaussian', 'Velocity Random walk'});




function h = plotPoses (poseCell)
    pos_arr = [];
    for ii = 1 : length(poseCell)
        pos_arr = [pos_arr, poseCell{ii}(:, 4)];
    end
    hold on;
    h = plot3(pos_arr(1, :), pos_arr(2, :), pos_arr(3, :), '.-');
    hold off;
    view(3);
end



