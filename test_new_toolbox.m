clear all
close all
clc


% Straight road — Straight Road
% 
% Curved road — Curved Road
% 
% Parking lot — Parking Lot
% 
% Double lane change — Double Lane Change
% 
% Open surface — Open Surface
% 
% US city block — US City Block
% 
% US highway — US Highway
% 
% Virtual Mcity — Virtual Mcity
% 
% Large parking lot — Large Parking Lot
% 


Simulation3DSensor



sceneName = 'LargeParkingLot';
[sceneImage, sceneRef] = helperGetSceneImage(sceneName);

hScene = figure;
helperShowSceneImage(sceneImage, sceneRef)
title(sceneName)


