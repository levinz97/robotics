#include <stdlib.h>
// #include <unistd.h>
#include <random>
#include <Algo/ann.h>
#include <Kin/kin.h>
#include <Kin/proxy.h>
#include <Gui/opengl.h>
#include <Kin/frame.h>
#include <Kin/viewer.h>
#include <Kin/F_collisions.h>

struct RRT{
private:
  ANN ann;      //ann stores all points added to the tree in ann.X
  uintA parent; //for each point we also store the index of the parent node
  double stepsize;
  uint nearest;

public:
  rai::Mesh lines;

public:
  RRT(const arr& q0, double _stepsize){
    ann.append(q0);      //append q as the root of the tree
    parent.append(0);    //q has itself as parent
    stepsize = _stepsize;
  }

  double getProposalTowards(arr& q){
    //find NN
    nearest=ann.getNN(q);

    //compute little step
    arr d = q - ann.X[nearest]; //difference vector between q and nearest neighbor
    double dist = length(d);
    q = ann.X[nearest] + stepsize/dist * d;
    return dist;
  }

  void add(const arr& q){
    ann.append(q);
    parent.append(nearest);
  }

  void addLineDraw(const arr& q, rai::Configuration& K){
    //We can't draw the edge in the 7-dim joint space!
    //But we can draw a projected edge in 3D endeffector position space:
    arr y_from,y_to;
    arr line;
    K.setJointState(ann.X[nearest]);  K.kinematicsPos(y_from, NoArr, K.getFrame("peg"));
    K.setJointState(q             );  K.kinematicsPos(y_to  , NoArr, K.getFrame("peg"));
    lines.V.append(y_from);
    lines.V.append(y_to);
    lines.V.reshape(lines.V.N/3, 3);
    lines.T.append({lines.V.d0-2, lines.V.d0-1});
    lines.T.reshape(lines.T.N/2, 2);
  }

  //some access routines
  uint getNearest(){ return nearest; }
  uint getParent(uint i){ return parent(i); }
  uint getNumberNodes(){ return ann.X.d0; }
  arr getNode(uint i){ return ann.X[i]; }
  void getRandomNode(arr& q){ q = ann.X[rnd(ann.X.d0)]; }
};


void RTTplan(){
  rai::Configuration C("pegInAHole.g");

  arr qT = {0.945499, 0.431195, -1.97155, 0.623969, 2.22355, -0.665206, -1.48356};
  arr q0, y_col, q;
  q0 = C.getJointState();
  // swap the q0 and qT according to 1.a
  q0 = qT.copy();
  qT = C.getJointState();

  q=q0;

  cout <<"final posture (hit ENTER in the OpenGL window to continue!!)" <<endl;
  C.setJointState(qT);
  C.watch(true);
  
  cout <<"initial posture (hit ENTER in the OpenGL window to continue!!)" <<endl;
  C.setJointState(q0);
  C.watch(true);
  
  C.watch(false);

  double stepsize = .1;
  RRT rrt0(q0, stepsize);

  rai::Frame *f = C.addFrame("lines0");
  f->setConvexMesh({});  // Add an empty mesh for line drawing to the world
  f->setContact(0);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.,1.);
  uint i;
  for(i=0;i<10000;i++){
    // if random_num < beta, set the q_goal(qT) to q_target
    if (dis(gen) < 0.5)
    {
      q = qT;
    } 
    else 
    {
      // else let rrt0 grow in random direction
      rndUniform(q,-RAI_2PI,RAI_2PI,false);
    }
    
    // compute q_new
    rrt0.getProposalTowards(q);
    C.setJointState(q);
    
    // check if q is collision free
    C.stepSwift();
    Value col = F_AccumulatedCollisions()
                .eval(C.frames);
    // C.kinematicsProxyCost(y_col, NoArr);
    if(col.y(0)<=1e-10){
      rrt0.add(q);
      rrt0.addLineDraw(q,C);
    }

    
    //some output
    if(!(i%100)){
      C["lines0"]->shape->mesh() = rrt0.lines;  // updates mesh lines0 with lines from rrt
      C.setJointState(q);
      C.gl()->recopyMeshes(C);
      C.watch(true);
      cout <<"\rRRT samples=" << i <<" tree sizes = " <<rrt0.getNumberNodes() << std::flush;
    }

      if (length(q-qT) < stepsize)
    {
      cout << "\nfind the solution !" << endl;
      break;
    }
    
  }
  C.watch(true);

  rrt0.lines.clear();

  uint nearest = rrt0.getNearest();
  arr  node = rrt0.getNode(nearest);
  uint node_idx = rrt0.getParent(nearest);

  f->setConvexMesh({});
  f->setContact(0);
  while (length(node - q0) >= stepsize)
  {
    node_idx = rrt0.getParent(node_idx);
    node = rrt0.getNode(node_idx);
    rrt0.addLineDraw(node, C);

    C["lines0"]->shape->mesh() = rrt0.lines;
    C.setJointState(node);
    C.gl()->recopyMeshes(C);
    C.watch(true);
    
  }
  
  
  C.watch(true);
}

///------------------------------------------------bidirectional_rrt-----------------------------------------------------
void RTTplan(bool use_bi_directional){
  rai::Configuration C("pegInAHole.g");

  arr qT = {0.945499, 0.431195, -1.97155, 0.623969, 2.22355, -0.665206, -1.48356};
  arr q0, y_col, q, q_; // q_ for backwards
  q0 = C.getJointState();

  q = q0;
  q_ = qT;

  cout <<"final posture (hit ENTER in the OpenGL window to continue!!)" <<endl;
  C.setJointState(qT);
  C.watch(true);
  
  cout <<"initial posture (hit ENTER in the OpenGL window to continue!!)" <<endl;
  C.setJointState(q0);
  C.watch(true);
  
  C.watch(false);

  double stepsize = .1;
  // implement bidirectional rrt
  RRT rrt_forward(q0, stepsize);
  RRT rrt_backward(qT, stepsize);

  rai::Frame *f = C.addFrame("lines0");
  f->setConvexMesh({});  // Add an empty mesh for line drawing to the world
  f->setContact(0);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.,1.);
  uint i;
  for(i=0;i<10000;i++){
    // if random_num < beta, set the q_goal(qT) to q_target
    if (dis(gen) < 0.5)
    {
      q = qT;
      q_ = q0;
    } 
    else 
    {
      // else let rrt0 grow in random direction
      rndUniform(q,-RAI_2PI,RAI_2PI,false);
      rndUniform(q_,-RAI_2PI,RAI_2PI,false);
    }
    // ----------------------------------------------------- forward part ---------------------------------------------------------------------------//
    // compute q_new
    rrt_forward.getProposalTowards(q);
    C.setJointState(q);
    
    // check if q is collision free
    C.stepSwift();
    Value col = F_AccumulatedCollisions()
                .eval(C.frames);
    // C.kinematicsProxyCost(y_col, NoArr);
    if(col.y(0)<=1e-10){
      rrt_forward.add(q);
      rrt_forward.addLineDraw(q,C);
    }
    
    //some output
    if(!(i%100)){
      C["lines0"]->shape->mesh() = rrt_forward.lines;  // updates mesh lines0 with lines from rrt
      C.setJointState(q);
      C.gl()->recopyMeshes(C);
      C.watch(true);
      cout <<"\rRRT_forward samples=" << i <<" tree sizes = " <<rrt_forward.getNumberNodes() << std::flush;
    }
    // ---------------------------------------------------------backward part-------------------------------------
    rrt_backward.getProposalTowards(q_);
    C.setJointState(q_);
    C.stepSwift();
    col = F_AccumulatedCollisions().eval(C.frames);
    if(col.y(0) <= 1e-10)
    {
      rrt_backward.add(q_);
      rrt_backward.addLineDraw(q_,C);
    }

  //some output
    if(!(i%100))
    {
      C["lines0"]->shape->mesh() = rrt_backward.lines;  // updates mesh lines0 with lines from rrt
      C.setJointState(q_);
      C.gl()->recopyMeshes(C);
      C.watch(true);
      cout <<"\rRRT_backward samples=" << i <<" tree sizes = " <<rrt_backward.getNumberNodes() << std::flush;
    }

    auto q_forw = q_.copy();
    auto q_backw = q.copy();

    rrt_forward.getProposalTowards(q_forw);
    if (length(q_forw - q_) < stepsize)
    {
      cout << "\nfind the solution !" << endl;
      break;
    }
    rrt_backward.getProposalTowards(q_backw);
    if (length(q-q_backw) < stepsize)
    {
      cout << "\nfind the solution !" << endl;
      break;
    }
    
  }
  C.watch(true);

  rrt_forward.lines.clear();
  rrt_backward.lines.clear();

  // rrt_forward.getProposalTowards(q0);
  uint nearest = rrt_forward.getNearest();
  arr  node = rrt_forward.getNode(nearest);
  uint node_idx = rrt_forward.getParent(nearest);

  f->setConvexMesh({});
  f->setContact(0);
  while (length(node - q) >= stepsize)
  {
    node_idx = rrt_forward.getParent(node_idx);
    node = rrt_forward.getNode(node_idx);
    rrt_forward.addLineDraw(node, C);

    C["lines0"]->shape->mesh() = rrt_forward.lines;
    C.setJointState(node);
    C.gl()->recopyMeshes(C);
    C.watch(true);
    
  }


  nearest = rrt_backward.getNearest();
  node = rrt_backward.getNode(nearest);
  node_idx = rrt_backward.getParent(nearest);

  f->setConvexMesh({});
  f->setContact(0);
  while (length(node - qT) >= stepsize)
  {
    node_idx = rrt_backward.getParent(node_idx);
    node = rrt_backward.getNode(node_idx);
    rrt_backward.addLineDraw(node, C);

    C["lines0"]->shape->mesh() = rrt_backward.lines;
    C.setJointState(node);
    C.gl()->recopyMeshes(C);
    C.watch(true);
    
  }
  cout <<"\rRRT_forward tree sizes = " <<rrt_forward.getNumberNodes() << std::endl;
  cout <<"\rRRT_backward tree sizes = " <<rrt_backward.getNumberNodes() << std::endl;
  
  C.watch(true);
}








int main(int argc,char **argv){
  rai::initCmdLine(argc,argv);

  // RTTplan();
  RTTplan(true);

  return 0;
}
