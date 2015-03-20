/*******************************************************************************
* Structured Edge Detection Toolbox      Version 3.01
* Code written by Piotr Dollar and Larry Zitnick, 2014.
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/
#include "mex.h"
#include "math.h"
#include <algorithm>
#include <vector>
#include <map>
using namespace std;
#define PI 3.14159265f
int clamp( int v, int a, int b ) { return v<a?a:v>b?b:v; }

// trivial array class encapsulating pointer arrays
template <class T> class Array
{
public:
  Array() { _h=_w=0; _x=0; _free=0; }
  virtual ~Array() { clear(); }
  void clear() { if(_free) delete [] _x; _h=_w=0; _x=0; _free=0; }
  void init(int h, int w) { clear(); _h=h; _w=w; _x=new T[h*w](); _free=1; }
  T& val(size_t c, size_t r) { return _x[c*_h+r]; }
  int _h, _w; T *_x; bool _free;
};

// convenient typedefs
typedef vector<float> vectorf;
typedef vector<int> vectori;
typedef Array<float> arrayf;
typedef Array<int> arrayi;

// bounding box data structures and routines
typedef struct { int c, r, w, h; float s; } Box;
typedef struct { int top, bottom, left, right;} myBox;
typedef vector<Box> Boxes;
bool boxesCompare( const Box &a, const Box &b ) { return a.s<b.s; }
float boxesOverlap( Box &a, Box &b );
void boxesNms( Boxes &boxes, float thr, int maxBoxes );

struct pairs {
  pairs(int segId, float weight) : segId(segId), weight(weight){}
  int segId;
  float weight;
};
bool pairsCompare( const pairs &a, const pairs &b ) { return a.weight<b.weight; }

//disjoint set class
typedef struct {
  int rank;
  int p;
  myBox boundingBox;
  float weight;
} uni_elt;

class universe {
public:
  universe(int elements, vector<myBox> _segmentBoxes);
  ~universe();
  float getMax();
  Box getBestBoundingBox();
  void updateWeights(int segId, float weight);
  void updateMax(int segId);
  int find(int x);
  void join(int x, int y);
  int computeArea(myBox boundingBox);
  myBox updateBoundingBox(myBox a, myBox b);
  void resetSegment(int x);
  void resetMax();
  int num_sets() const {return num; }

private:
  uni_elt *elts;
  int num;
  float maximumScore;
  vector<myBox> segmentBoxes;
  myBox bestBoundingBox;
};

universe::universe(int elements, vector<myBox> _segmentBoxes) {
  elts = new uni_elt[elements];
  num = elements;
  segmentBoxes = _segmentBoxes;
}

universe::~universe() {
  delete [] elts;
}

void universe::resetSegment(int x){
  elts[x].rank = 0;
  elts[x].p = x;
  elts[x].boundingBox = segmentBoxes[x];
  elts[x].weight = 0;
}

void universe::resetMax(){
  maximumScore = 0;
}

float universe::getMax() {
  return maximumScore;
}

Box universe::getBestBoundingBox(){
  Box x;
  x.c = bestBoundingBox.left;
  x.r = bestBoundingBox.top;
  x.w = bestBoundingBox.right - bestBoundingBox.left;
  x.h = bestBoundingBox.bottom - bestBoundingBox.top;
  x.s = getMax();
  return x;
}

void universe::updateMax(int segId) {
  int area = computeArea(elts[segId].boundingBox);
  float score = area*elts[segId].weight;
  if (maximumScore<score){
    maximumScore = score;
    bestBoundingBox = elts[segId].boundingBox;
  }
}

void universe::updateWeights(int segId, float weight) {
  elts[segId].weight = weight;
}

int universe::computeArea(myBox box) {
  return ((box.bottom - box.top) * (box.right - box.left));
}

myBox universe::updateBoundingBox(myBox a, myBox b) {
  myBox resultBox;
  resultBox.top = min(a.top,b.top);
  resultBox.bottom = max(a.bottom,b.bottom);
  resultBox.left = min(a.left,b.left);
  resultBox.right = max(a.right,b.right);
  return resultBox;
}

int universe::find(int x) {
  int y = x;
  while (y != elts[y].p)
    y = elts[y].p;
  elts[x].p = y;
  return y;
}

void universe::join(int x, int y) {

  if (elts[x].rank > elts[y].rank) {
    elts[y].p = x;
    elts[x].boundingBox = updateBoundingBox(elts[x].boundingBox,elts[y].boundingBox);
    elts[x].weight = min(elts[x].weight,elts[y].weight);
    updateMax(x);
  }
  else {
    elts[x].p = y;
    elts[y].boundingBox = updateBoundingBox(elts[x].boundingBox,elts[y].boundingBox);
    elts[y].weight = min(elts[x].weight,elts[y].weight);
    updateMax(y);

    if (elts[x].rank == elts[y].rank)
      elts[y].rank++;
  }
}

// main class for generating edge boxes
class SegmentBoxGenerator
{
public:
  // method parameters (must be manually set)
  float _alpha, _beta, _minScore; int _maxBoxes;
  float _edgeMinMag, _edgeMergeThr, _clusterMinMag;
  float _maxAspectRatio, _minBoxArea, _gamma, _kappa;

  // main external routine (set parameters first)
  void generate( Boxes &boxes, arrayf &I, arrayf &edges );

private:
  // edge segment information (see clusterEdges)
  int h, w;                         // image dimensions
  int _segCnt;                      // total segment count
  arrayi _segIds;                   // segment ids (-1/0 means no segment)
  vectorf _segMag;                  // segment edge magnitude sums
  arrayf E1;
  vector<vectorf> _segAff;          // segment affinities
  vector<vectori> _segAffIdx;       // segment neighbors
  vector<myBox> _segmentBoxes;

  vectori _segR, _segC;             // segment lower-right pixel
  // data structures for efficiency (see prepDataStructs)
  arrayf _magIImg; arrayi _hIdxImg, _vIdxImg;
  vector<vectori> _hIdxs, _vIdxs; vectorf _scaleNorm;
  float _scStep, _arStep, _rcStepRatio;

  // data structures for efficiency (see scoreBox)
  arrayf _sWts, _sDist; arrayi _sDone, _sMap, _sIds, _sBlacked; int _sId;

  universe *u;

  // helper routines
  void prepDataStructs( arrayf &I );
  void scoreAllBoxes( Boxes &boxes );
  void scoreBox( Box &box );
  void refineBox( Box &box );
  void createSegments( arrayf &I, arrayf &edges );
  void createSegmentBoxes( arrayf &I );
};

////////////////////////////////////////////////////////////////////////////////

void SegmentBoxGenerator::generate( Boxes &boxes, arrayf &I, arrayf &edges )
{
  createSegments(I,edges);
  createSegmentBoxes(I);
  prepDataStructs(I); 
  scoreAllBoxes(boxes);
}

void SegmentBoxGenerator::createSegments( arrayf &I, arrayf &edges )
{
    int c, r, j; h=I._h; w=I._w;
    vectori _R, _C;

    // with image borders

    // w=w+2;
    // h=h+2;

    // _segIds.init(h,w);
    // for( c=0; c<w; c++ ) for( r=0; r<h; r++ ) {
    //     if( c==0 || r==0 || c==w-1 || r==h-1)
    //         _segIds.val(c,r)=-1;
    //     else {
    //         _segIds.val(c,r)=int(I.val(c-1,r-1));
    //         if (I.val(c-1,r-1) > _segCnt)
    //             _segCnt = int(I.val(c-1,r-1));
    //     }
    // }
    // _segCnt++; // because we need to count the 0 id segments
    
    // without image borders
    
    _segIds.init(h,w);
       for( c=0; c<w; c++ ) for( r=0; r<h; r++ ) {
        _segIds.val(c,r)=int(I.val(c,r));
        if (I.val(c,r) > _segCnt)
            _segCnt = int(I.val(c,r));
    }
    _segCnt++; // because we need to count the 0 id segments
     
    // segments affinities
    _segAff.resize(_segCnt); _segAffIdx.resize(_segCnt);
    
    float normalize = float(sqrt(3.0)/256);
    
    for(int i=0; i<edges._h; i++){
      int a = (int)edges.val(0,i);
      int b = (int)edges.val(1,i);
      float w = (edges.val(2,i))*normalize; //kanonikopoiisi sto 1 
      _segAff[a].push_back(w); _segAffIdx[a].push_back(b);
      _segAff[b].push_back(w); _segAffIdx[b].push_back(a);
    }

    // compute _segC and _segR
    _segC.resize(_segCnt); _segR.resize(_segCnt);
    for( c=0; c<w; c++ ) for( r=0; r<h; r++ ) {
      if( (j=_segIds.val(c,r))>-1) {
        _segC[j]=c; _segR[j]=r;
      }
    }
}

void SegmentBoxGenerator::prepDataStructs( arrayf &I )
{
  int c, r, i;

  // initialize step sizes
  _scStep=sqrt(1/_alpha);
  _arStep=(1+_alpha)/(2*_alpha);
  _rcStepRatio=(1-_alpha)/(1+_alpha);

  // create _scaleNorm
  _scaleNorm.resize(10000);
  for( i=0; i<10000; i++ )
    _scaleNorm[i]=pow(1.f/i,_kappa);

  // create remaining data structures
  // den elegxei an pernaei to idio segment apo to orio. Logika einai
  // afou dinei akoma mikroteri timi gia to score  
  _hIdxs.resize(h); _hIdxImg.init(h,w);
  for( r=0; r<h; r++ ) {
    int s=0, s1; _hIdxs[r].push_back(s);
    for( c=0; c<w; c++ ) {
      s1 = _segIds.val(c,r);
      if( s1!=s ) { s=s1; _hIdxs[r].push_back(s); }
      _hIdxImg.val(c,r) = int(_hIdxs[r].size())-1;
    }
  }
  _vIdxs.resize(w); _vIdxImg.init(h,w);
  for( c=0; c<w; c++ ) {
    int s=0; _vIdxs[c].push_back(s);
    for( r=0; r<h; r++ ) {
      int s1 = _segIds.val(c,r);
      if( s1!=s ) { s=s1; _vIdxs[c].push_back(s); }
      _vIdxImg.val(c,r) = int(_vIdxs[c].size())-1;
    }
  }

  // initialize scoreBox() data structures
  int n=_segCnt+1; 
  _sWts.init(n,1); _sDist.init(n,1);
  _sBlacked.init(n,1); _sDone.init(n,1); 
  _sMap.init(n,1); _sIds.init(n,1);
  for( i=0; i<n; i++ ) {
    _sDone.val(0,i)=-1; 
    _sBlacked.val(0,i)=-1;
  }
  _sId=0;
}

void SegmentBoxGenerator::createSegmentBoxes( arrayf &I )
{
  int c, r, j;

  _segmentBoxes.resize(_segCnt);
  vectori alreadyDone;
  alreadyDone.resize(_segCnt,0);

  for( c=0; c<w; c++ ) for( r=0; r<h; r++ ){
    if ((j=_segIds.val(c,r))>=0){  
      if (alreadyDone[j]==0){
        _segmentBoxes[j].top = r;
        _segmentBoxes[j].bottom = r+1;
        _segmentBoxes[j].left = c;
        _segmentBoxes[j].right = c+1;
        alreadyDone[j] = 1;
      }
      else {
        if ( c < _segmentBoxes[j].left ) {
          _segmentBoxes[j].left = c;
        }
        if ( c > _segmentBoxes[j].right ) {
          _segmentBoxes[j].right = c+1;
        }
        if ( r < _segmentBoxes[j].top ) {
          _segmentBoxes[j].top = r;
        }
        if ( r > _segmentBoxes[j].bottom ) {
          _segmentBoxes[j].bottom = r+1;
        }
      }
    }
  }
  u = new universe(_segCnt, _segmentBoxes);
}

void SegmentBoxGenerator::scoreBox( Box &box )
{
  int i, j, k, q, bh, bw, r0, c0, r1, c1;
  float *sWts=_sWts._x, *sDist=_sDist._x; 
  int sId=_sId++;
  int *sDone=_sDone._x, *sMap=_sMap._x, *sIds=_sIds._x;
  int *sBlacked=_sBlacked._x;

  // add edge count inside box
  r1=clamp(box.r+box.h,0,h-1); r0=box.r=clamp(box.r,0,h-1);
  c1=clamp(box.c+box.w,0,w-1); c0=box.c=clamp(box.c,0,w-1);
  bh=box.h=r1-box.r; bw=box.w=c1-box.c; 
  // bw/=2;bh/=2; 
  
  float norm = pow(1.f/(bw*bh),_kappa);
  
  // find interesecting segments along four boundaries
  
  int cs, ce, rs, re, n=0;
  cs=_hIdxImg.val(c0,r0); ce=_hIdxImg.val(c1,r0); // top
  for( i=cs; i<=ce; i++ ) if( (j=_hIdxs[r0][i])>=0 && sDone[j]!=sId ) {
    sIds[n]=j; sWts[n]=1; sDist[n]=0; sDone[j]=sId; sMap[j]=n++;
  }
  cs=_hIdxImg.val(c0,r1); ce=_hIdxImg.val(c1,r1); // bottom
  for( i=cs; i<=ce; i++ ) if( (j=_hIdxs[r1][i])>=0 && sDone[j]!=sId ) {
    sIds[n]=j; sWts[n]=1; sDist[n]=0; sDone[j]=sId; sMap[j]=n++;
  }
  rs=_vIdxImg.val(c0,r0); re=_vIdxImg.val(c0,r1); // left
  for( i=rs; i<=re; i++ ) if( (j=_vIdxs[c0][i])>=0 && sDone[j]!=sId ) {
    sIds[n]=j; sWts[n]=1; sDist[n]=0; sDone[j]=sId; sMap[j]=n++;
  }
  rs=_vIdxImg.val(c1,r0); re=_vIdxImg.val(c1,r1); // right
  for( i=rs; i<=re; i++ ) if( (j=_vIdxs[c1][i])>=0 && sDone[j]!=sId ) {
    sIds[n]=j; sWts[n]=1; sDist[n]=0; sDone[j]=sId; sMap[j]=n++;
  }
  
  // follow connected paths and set weights accordingly (w=0 means remove)
  vector<pairs> segmentsInside;

  for( i=0; i<n; i++ ) {
    float w=sDist[i]; j=sIds[i];
    for( k=0; k<int(_segAffIdx[j].size()); k++ ) {
      q=_segAffIdx[j][k];
      float wq=max(w,_segAff[j][k]); // because big segAff means big difference
      // wq = 1;
      if( sDone[q]==sId ) {
        if( wq<sDist[sMap[q]] ) { sDist[sMap[q]]=wq; i=min(i,sMap[q]-1); }
      } 
      else if(_segC[q]>=c0 && _segC[q]<=c1 && _segR[q]>=r0 && _segR[q]<=r1) {
        sIds[n]=q; sDist[n]=wq; sWts[n]=0; sDone[q]=sId; sMap[q]=n++;
        u->resetSegment(q);
      }
    }
  }

  u->resetMax();

  for ( i=0; i<n; i++ ) {
    if (sDist[i]>0) {
      segmentsInside.push_back(pairs(sIds[i],sDist[i]));
      u->updateWeights(sIds[i],sDist[i]);
      u->updateMax(sIds[i]);
    }
  }

  sort(segmentsInside.rbegin(),segmentsInside.rend(),pairsCompare);

  for ( i=0; i<segmentsInside.size(); i++ ) {
    j = segmentsInside[i].segId;
    sBlacked[j] = sId;
    for ( k=0; k<int(_segAffIdx[j].size()); k++) {
      q=_segAffIdx[j][k];
      if (sBlacked[q] == sId) {
        int a = u->find(j);
        int b = u->find(q);
        if (a!=b) {
          u->join(a,b);
        }
      }
    }
  }
  // mexPrintf("%f\n", u->getMax());

  // float v = u->getMax();

  // v*=norm;
  // if(v<_minScore) v=0; 
  // box.s=v;

  box = u->getBestBoundingBox();
  box.s *= norm;
}

void SegmentBoxGenerator::refineBox( Box &box )
{
  int rStep = int(box.h*_rcStepRatio);
  int cStep = int(box.w*_rcStepRatio);
  while( 1 ) {
    // prepare for iteration
    rStep/=2; cStep/=2; if( rStep<=2 && cStep<=2 ) break;
    rStep=max(1,rStep); cStep=max(1,cStep); Box B;
    // search over r start
    B=box; B.r=box.r-rStep; B.h=B.h+rStep; scoreBox(B);
    if(B.s<=box.s) { B=box; B.r=box.r+rStep; B.h=B.h-rStep; scoreBox(B); }
    if(B.s>box.s) box=B;
    // search over r end
    B=box; B.h=B.h+rStep; scoreBox(B);
    if(B.s<=box.s) { B=box; B.h=B.h-rStep; scoreBox(B); }
    if(B.s>box.s) box=B;
    // search over c start
    B=box; B.c=box.c-cStep; B.w=B.w+cStep; scoreBox(B);
    if(B.s<=box.s) { B=box; B.c=box.c+cStep; B.w=B.w-cStep; scoreBox(B); }
    if(B.s>box.s) box=B;
    // search over c end
    B=box; B.w=B.w+cStep; scoreBox(B);
    if(B.s<=box.s) { B=box; B.w=B.w-cStep; scoreBox(B); }
    if(B.s>box.s) box=B;
  }
}

void SegmentBoxGenerator::scoreAllBoxes( Boxes &boxes )
{
  // get list of all boxes roughly distributed in grid
  boxes.resize(0); int arRad, scNum; float minSize=sqrt(_minBoxArea);
  arRad = int(log(_maxAspectRatio)/log(_arStep*_arStep));
  scNum = int(ceil(log(max(w,h)/minSize)/log(_scStep)));
  for( int s=0; s<scNum; s++ ) {
    int a, r, c, bh, bw, kr, kc, bId=-1; float ar, sc;
    for( a=0; a<2*arRad+1; a++ ) {
      ar=pow(_arStep,float(a-arRad)); sc=minSize*pow(_scStep,float(s));
      bh=int(sc/ar); kr=max(2,int(bh*_rcStepRatio));
      bw=int(sc*ar); kc=max(2,int(bw*_rcStepRatio));
      for( c=0; c<w-bw+kc; c+=kc ) for( r=0; r<h-bh+kr; r+=kr ) {
        Box b; b.r=r; b.c=c; b.h=bh; b.w=bw; boxes.push_back(b);
      }
    }
  }

  // boxes.resize(0);
  // Box b; 
  // b.r=135; b.c=746; b.h = 127; b.w= 151;
  // boxes.push_back(b);
  // b.r=235; b.c=652; b.h = 110; b.w= 120;
  // boxes.push_back(b);

  // score all boxes, refine top candidates, perform nms
  int i, k=0, m = int(boxes.size());
  for( i=0; i<m; i++ ) {
    scoreBox(boxes[i]);
    // mexPrintf("\nnew one\n");
    if( !boxes[i].s ) continue; k++;
    refineBox(boxes[i]);
  }
  sort(boxes.rbegin(),boxes.rend(),boxesCompare);
  boxes.resize(k); boxesNms(boxes,_beta,_maxBoxes);
}

float boxesOverlap( Box &a, Box &b ) {
  float areai, areaj, areaij;
  int r0, r1, c0, c1, r1i, c1i, r1j, c1j;
  r1i=a.r+a.h; c1i=a.c+a.w; if( a.r>=r1i || a.c>=c1i ) return 0;
  r1j=b.r+b.h; c1j=b.c+b.w; if( a.r>=r1j || a.c>=c1j ) return 0;
  areai = (float) a.w*a.h; r0=max(a.r,b.r); r1=min(r1i,r1j);
  areaj = (float) b.w*b.h; c0=max(a.c,b.c); c1=min(c1i,c1j);
  areaij = (float) max(0,r1-r0)*max(0,c1-c0);
  return areaij / (areai + areaj - areaij);
}

void boxesNms( Boxes &boxes, float thr, int maxBoxes )
{
  sort(boxes.rbegin(),boxes.rend(),boxesCompare);
  if( thr>.99 ) return; const int nBin=10000;
  const float step=1/thr, lstep=log(step);
  vector<Boxes> kept; kept.resize(nBin+1);
  int i=0, j, k, n=(int) boxes.size(), m=0, b;
  while( i<n && m<maxBoxes ) {
    b = boxes[i].w*boxes[i].h; bool keep=1;
    b = clamp(int(ceil(log(float(b))/lstep)),1,nBin-1);
    for( j=b-1; j<=b+1; j++ )
      for( k=0; k<kept[j].size(); k++ ) if( keep )
        keep = boxesOverlap( boxes[i], kept[j][k] ) <= thr;
    if(keep) { kept[b].push_back(boxes[i]); m++; } i++;
  }
  boxes.resize(m); i=0;
  for( j=0; j<nBin; j++ )
    for( k=0; k<kept[j].size(); k++ )
      boxes[i++]=kept[j][k];
  sort(boxes.rbegin(),boxes.rend(),boxesCompare);
}

////////////////////////////////////////////////////////////////////////////////

// Matlab entry point: bbs = mex( E, O, prm1, prm2, ... )
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  // check and get inputs
  if(nr != 13) mexErrMsgTxt("Thirteen inputs required.");
  if(nl > 2) mexErrMsgTxt("At most two outputs expected.");
  if(mxGetClassID(pr[0])!=mxSINGLE_CLASS) mexErrMsgTxt("I must be a float*");
  arrayf I; I._x = (float*) mxGetData(pr[0]);
  arrayf edges; edges._x = (float*)mxGetData(pr[1]);
  int h = (int) mxGetM(pr[0]); I._h=h;
  int w = (int) mxGetN(pr[0]); I._w=w;
  edges._h = (int) mxGetM(pr[1]);
  edges._w = (int) mxGetN(pr[1]);

  // setup and run SegmentBoxGenerator
  SegmentBoxGenerator segmentBoxGen; Boxes boxes;
  
  segmentBoxGen._alpha = float(mxGetScalar(pr[2]));
  segmentBoxGen._beta = float(mxGetScalar(pr[3]));
  segmentBoxGen._minScore = float(mxGetScalar(pr[4]));
  segmentBoxGen._maxBoxes = int(mxGetScalar(pr[5]));
  segmentBoxGen._edgeMinMag = float(mxGetScalar(pr[6]));
  segmentBoxGen._edgeMergeThr = float(mxGetScalar(pr[7]));
  segmentBoxGen._clusterMinMag = float(mxGetScalar(pr[8]));
  segmentBoxGen._maxAspectRatio = float(mxGetScalar(pr[9]));
  segmentBoxGen._minBoxArea = float(mxGetScalar(pr[10]));
  segmentBoxGen._gamma = float(mxGetScalar(pr[11]));
  segmentBoxGen._kappa = float(mxGetScalar(pr[12]));
   
  segmentBoxGen.generate( boxes, I, edges );

  // create output bbs and output to Matlab
  int n = (int) boxes.size();
  pl[0] = mxCreateNumericMatrix(n,5,mxSINGLE_CLASS,mxREAL);
  float *bbs = (float*) mxGetData(pl[0]);
  for(int i=0; i<n; i++) {
    bbs[ i + 0*n ] = (float) boxes[i].c+1;
    bbs[ i + 1*n ] = (float) boxes[i].r+1;
    bbs[ i + 2*n ] = (float) boxes[i].w;
    bbs[ i + 3*n ] = (float) boxes[i].h;
    bbs[ i + 4*n ] = boxes[i].s;
 }
}
