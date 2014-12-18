function bbs = dokimi(dokimastiko)

alpha = 1;
beta= 1;
minScore= 1;
maxBoxes= 1;
segmentMinMag= 1;
segmentMergeThr= 1;
clusterMinMag= 1;
maxAspectRatio= 1;
minBoxArea= 1;
gamma= 1;
kappa= 1;

bbs=segmentBoxesMex(dokimastiko,alpha,beta,minScore,maxBoxes,...
  segmentMinMag,segmentMergeThr,clusterMinMag,...
  maxAspectRatio,minBoxArea,gamma,kappa);

end