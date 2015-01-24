%% set up opts for segmentBoxes (see segmentBoxes.m)
opts = segmentBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

% %% detect Edge Box bounding box proposals (see segmentBoxes.m)
% 
% tic, bbs=segmentBoxes('peppers.txt','edges.txt','peppers',opts); toc % 'peppers' is the id of the image
% 
% %% show evaluation results (using pre-defined or interactive boxes)
% gt=[122 248 92 65; 193 82 71 53; 410 237 101 81; 204 160 114 95; ...
%   9 185 86 90; 389 93 120 117; 253 103 107 57; 81 140 91 63; 77 100 107 63; 309 164 191 133];
%  gt = bbs(1:16,1:4);
% % gt = bbs(1,1:4);
% % gt = bbs(50:70,1:4);
% if(0), gt='Please select an object box.'; disp(gt); figure(1); imshow(I);
%   title(gt); [~,gt]=imRectRot('rotate',0); gt=gt.getPos(); end
% gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
% original = imread('peppers.png');
% figure(1); bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
% title('green=matched gt  red=missed gt  dashed-green=matched detect');
% 
% %% run and evaluate on entire dataset (see boxesData.m and boxesEval.m)
% % if(~exist('boxes/VOCdevkit/','dir')), return; end
% % split='val'; data=boxesData('split',split);
% % nm='segmentBoxes70'; opts.name=['boxes/' nm '-' split '.mat'];
% % segmentBoxes(data.imgs,data.edges,data.ids,opts); opts.name=[];
% % boxesEval('data',data,'names',nm,'thrs',.7,'show',2);
% % boxesEval('data',data,'names',nm,'thrs',.5:.05:1,'cnts',1000,'show',3);

%% Pascal images test
imgdir = 'boxes\VOCdevkit\VOC2007\JPEGImages\';
txtdir = 'boxes\VOCdevkit\Segments\Segments\';
edgedir = 'boxes\VOCdevkit\Edges\Edges\';

%%% to change
testImage = '000002';
figure(1);

bbs=segmentBoxes([txtdir '' testImage '.txt'],[edgedir '' testImage '.txt'],testImage,opts);

split='val'; data=boxesData('split',split);
gt = data.gt;
n = data.n;
gt = gt(1:n);

%%% to change
gt = cell2mat(gt(2));

gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
original = imread([imgdir '' testImage '.jpg']);
subplot(1,2,1);
bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
title(['green=matched gt' char(10) 'red=missed gt' char(10) 'dashed-green=matched detect']);

%%% top scoring boxes
gt = bbs(1:10,1:4);

gt(:,5)=0; [gtRes,dtRes]=bbGt('evalRes',gt,double(bbs),.7);
original = imread([imgdir '' testImage '.jpg']);
subplot(1,2,2);
bbGt('showRes',original,gtRes,dtRes(dtRes(:,6)==1,:));
title('10 top scoring boxes');