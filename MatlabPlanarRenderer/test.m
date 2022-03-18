clear all
close all
clc


    %
    % USAGE: 
    tex = imread('licensePlates.jpg');
    mPR = MatlabPlanarRenderer(tex);
    T = [1 0 0 -150; 0 1 0 -150; 0 0 1 1100.5; 0 0 0 1];
    R = rotx(30)*roty(40)*rotz(25);
    T(1:3,1:3) = R;
    mPR.setTransformation(T);
    [img] = mPR.render();
    imshow(img)
    %    
    %   Point correspondences can be obtained from: mPR.getCorrespondence()
    %
    %   Also check the usage of:  setIntrinsics, toggleVerbosity,
    %   setDensity, getCorrespondence
    %
    %   Point correspondence format is described in the code comments

