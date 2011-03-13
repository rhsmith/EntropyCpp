clc;

load('entrop.mat');
load('pressureTotal.mat');
load('volumetotal.mat');
load('centrop.txt');
load('cpressure.txt');
load('cvolume.txt');

entropdiff = abs(entrop - centrop);
volumediff = abs(volumetotal - cvolume);
pressurediff = abs(pressuretotal - cpressure);

max(max(entropdiff))
max(max(volumediff))
max(max(pressurediff))

cvolume(:,1:4)
volumetotal(:,1:4)
