clc;clear;

%% set up opts for segmentBoxes (see segmentBoxes.m)
opts = segmentBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% detect Edge Box bounding box proposals (see segmentBoxes.m)
Original = imread('peppers.png');
I = imread('peppersSeg.ppm');
I = rgb2gray(I);
I = single(I);
I = I/255; %gia na kanonikopoithei
tic, bbs=segmentBoxes(I,opts); toc

%% show evaluation results (using pre-defined or interactive boxes)
gt=[122 248 92 65; 193 82 71 53; 410 237 101 81; 204 160 114 95; ...
  9 185 86 90; 389 93 120 117; 253 103 107 57; 81 140 91 63];
if(0), gt='Please select an object box.'; disp(gt); figure(1); imshow(I);
  title(gt); [~,gt]=imRectRot('rotate',0); gt=gt.getPos(); end
gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
figure(1); bbGt('showRes',Original,gtRes,dtRes(dtRes(:,6)==1,:));
title('green=matched gt  red=missed gt  dashed-green=matched detect');

%% run and evaluate on entire dataset (see boxesData.m and boxesEval.m)
%if(~exist('boxes/VOCdevkit/','dir')), return; end
%split='val'; data=boxesData('split',split);
%nm='segmentBoxes70'; opts.name=['boxes/' nm '-' split '.mat'];
%segmentBoxes(data.imgs,model,opts); opts.name=[];
%boxesEval('data',data,'names',nm,'thrs',.7,'show',2);
%boxesEval('data',data,'names',nm,'thrs',.5:.05:1,'cnts',1000,'show',3);
