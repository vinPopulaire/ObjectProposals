function bbs = dokimi(I)

I = rgb2gray(I);
I = single(I);
I = I/255; %gia na kanonikopoithei

alpha = .65;
beta= .75;
minScore= .01;
maxBoxes= 1e4;
segmentMinMag= .1;
segmentMergeThr= .5;
clusterMinMag= .5;
maxAspectRatio= 3;
minBoxArea= 1000;
gamma= 2;
kappa= 1.5;

bbs=segmentBoxesMex(I,alpha,beta,minScore,maxBoxes,...
  segmentMinMag,segmentMergeThr,clusterMinMag,...
  maxAspectRatio,minBoxArea,gamma,kappa);

end