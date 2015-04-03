%% load pre-trained edge detection model and set opts (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;
model.opts.multiscale=0; model.opts.sharpen=2; model.opts.nThreads=4;

%% set up opts for segmentBoxes (see segmentBoxes.m)
opts = segmentBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% test image
% tic, bbs=segmentBoxes('test\testout.txt','test\testedges.txt','test',opts); toc % 'peppers' is the id of the image
% figure(1);
% 
% gt = bbs(1:8,1:4);
% gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
% original = imread('test\test.pgm'); 
% subplot(1,2,1);
% bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
% title('green=matched gt  red=missed gt  dashed-green=matched detect');
% 
% original = imread('test\testseg.ppm'); 
% subplot(1,2,2);
% bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
% title('green=matched gt  red=missed gt  dashed-green=matched detect');



%% detect Edge Box bounding box proposals (see segmentBoxes.m)

% tic, bbs=segmentBoxes('peppers.txt','edges.txt','peppers',opts); toc % 'peppers' is the id of the image
% 
% %% show evaluation results (using pre-defined or interactive boxes)
% gt=[122 248 92 65; 193 82 71 53; 410 237 101 81; 204 160 114 95; ...
%   9 185 86 90; 389 93 120 117; 253 103 107 57; 81 140 91 63; 77 100 107 63; 309 164 191 133];
% % gt = bbs(1:16,1:4);
% if(0), gt='Please select an object box.'; disp(gt); figure(1); imshow(I);
%   title(gt); [~,gt]=imRectRot('rotate',0); gt=gt.getPos(); end
% gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
% original = imread('peppers.png');
% figure(1); bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
% title('green=matched gt  red=missed gt  dashed-green=matched detect');

%% run and evaluate on entire dataset (see boxesData.m and boxesEval.m)
% if(~exist('boxes/VOCdevkit/','dir')), return; end
% split='val'; data=boxesData('split',split);
% nm='segmentBoxes70'; opts.name=['boxes/' nm '-' split '.mat'];
% segmentBoxes(data.imgs,opts); opts.name=[];
% boxesEval('data',data,'names',nm,'thrs',.7,'show',2);
% boxesEval('data',data,'names',nm,'thrs',.5:.05:1,'cnts',1000,'show',3);

%% Pascal images test
imgdir = '..\ObjectProposals\boxes\VOCdevkit\VOC2007\JPEGImages\';

%% to change
testImage = '000003';
% figure(1);
I = imread([imgdir '' testImage '.jpg']);

tic, bbs=segmentBoxes(I,model,opts); toc;

f = figure(2);

split='val'; data=boxesData('split',split);
gt = data.gt;
n = data.n;
gt = gt(1:n);

%%% to change
gt = cell2mat(gt(3));

gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
subplot(1,3,1);
bbGt('showRes',I,gtRes,dtRes(dtRes(:,6)==1,:));
title(['green=matched gt' char(10) 'red=missed gt' char(10) 'dashed-green=matched detect']);

%%% top scoring boxes
gt = bbs(1:8,1:4);

gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
subplot(1,3,2);
bbGt('showRes',I,gtRes,dtRes(dtRes(:,6)==1,:));
title('8 top scoring boxes');

% gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
% original = imread([segmentdir '' testImage '.pgm']);
% subplot(1,3,3);
% bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
% title('8 top scoring boxes');

% %%% to change
% print -dpng 000011;

% figure(2);
% gt = bbs(1:2,1:4);
% gt(1,:) = bbs(1,1:4);
% gt(2,:) = bbs(6,1:4);
% gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
% original = imread([imgdir '' testImage '.jpg']);
% bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
% title('8 top scoring boxes');

fclose all;