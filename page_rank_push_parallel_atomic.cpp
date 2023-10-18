#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <atomic>
#define DEFAULT_STRATEGY "1"
#define DEFAULT_GRANULARITY "1"

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef double PageRankType;
#endif
static PageRankType *pr_curr;
static PageRankType *pr_next;
static double *timeTaken;
static double *barrier1_time;
static double *barrier2_time; 
static double* nextVertexToBeProcessed_time; 
std::atomic <PageRankType> *atomics_current;
std::atomic<PageRankType> *atomics_next;


typedef struct worker{
    uint threadID;
    Graph* graph_ptr;
    uint max_iters;
    uint total_num_vertices;
    CustomBarrier* barrier;
    uint startVertex;
    uint endVertex;
    double *time_buffer;
    PageRankType* pr_curr;
    PageRankType* pr_next;
    uint granularity;
    std::atomic <PageRankType> *atomics;
    int edgesProcessed;
    int verticesProcessed;
    double *barrier1_time;
    double *barrier2_time;
    double *nextVertexToBeProcessed_time;
}worker;
void Distributed_PageRank_Calculation(worker* worker)
{

  timer t;
  double time_taken=0.0;
  PageRankType* prr_curr_local= worker->pr_curr;
  PageRankType* prr_next_local= worker->pr_next;  
  uint start_vertex= worker->startVertex;
  uint end_vertex= worker->endVertex;
  uint max_iters= worker->max_iters;
  CustomBarrier *barrier_local= worker->barrier;
  Graph* g= worker->graph_ptr;
  std:: atomic <PageRankType>* atomics_local= worker->atomics;
  t.start();
  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = start_vertex; u <= end_vertex; u++) {
      uintE out_degree = g->vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g->vertices_[u].getOutNeighbor(i);

        //compare and swap here
        PageRankType fetch_old = atomics_next[v].load();

        PageRankType sum = fetch_old+(atomics_current[u].load() / (PageRankType) out_degree);

        while(!atomics_next[v].compare_exchange_weak(fetch_old, sum))
        {
          sum= fetch_old+(atomics_current[u].load() / (PageRankType) out_degree);
        }
        
      }
    }
    barrier_local->wait();
    
    for (uintV v = start_vertex; v <= end_vertex; v++) {
      atomics_next[v].store(PAGE_RANK(atomics_next[v].load()));

      // reset pr_curr for the next iteration
      atomics_current[v].store(atomics_next[v].load());
      atomics_next[v].store (0.0);
    }
    barrier_local->wait();

  }
  time_taken=t.stop();
  worker->time_buffer[worker->threadID]+=time_taken;
}

void pageRank_Init(Graph &g)
{
    uintV n = g.n_;
    pr_curr= new PageRankType[n];
    pr_next= new PageRankType[n];
    for(uintV i=0; i<n; i++)
    {
        pr_curr[i]=INIT_PAGE_RANK;
        pr_next[i]=0.0;        
    }
}
void pageRank_Destroy()
{
    delete[] pr_curr;
    delete[] pr_next;
}
void strategy1(uint &n_threads, Graph&g, uint &max_iterations)
{
  CustomBarrier* barrier= new CustomBarrier(n_threads);    
  //Initialzing array of doubles to keep track of time for each thread
  timeTaken= new double[n_threads]();
  barrier1_time= new double[n_threads]();
  barrier2_time= new double[n_threads]();
  nextVertexToBeProcessed_time= new double[n_threads]();
  // Initializing worker threads
  std::thread threads[n_threads];
  worker workers[n_threads];
//splitting the number of pages between each thread(evenly)
    
  std::vector<uint> startx(n_threads);
  std::vector<uint> endx(n_threads);
  
  uint min_pages_for_each_threads = g.n_ /  n_threads;
  uint excess_columns = g.n_ % n_threads;
  uint curr_page = 0;
  for (uint i = 0; i < n_threads; i++) {
    startx[i] = curr_page;
    if (excess_columns > 0) {
      endx[i] = curr_page + min_pages_for_each_threads;
      excess_columns--;
      } 
    else {
           endx[i] = curr_page + min_pages_for_each_threads - 1;

      }
    curr_page = endx[i]+1;
  } 
  pageRank_Init(g);
  //mutex for each vertex

  //std::mutex* mtxs= (std::mutex*) malloc(sizeof(std::mutex)* g.n_);
  atomics_current= new std::atomic<PageRankType>[g.n_];
  atomics_next= new std::atomic<PageRankType>[g.n_];

  for(int i=0;i< n_threads;i++)
  {
    workers[i].barrier= barrier;
    workers[i].startVertex= startx[i];
    workers[i].endVertex= endx[i];
    workers[i].graph_ptr= &g;
    workers[i].pr_curr= pr_curr;
    workers[i].pr_next= pr_next;
    workers[i].threadID=i;
    workers[i].time_buffer=timeTaken;
    workers[i].total_num_vertices=g.n_;
    workers[i].max_iters= max_iterations;
    workers[i].atomics= atomics_current;
    workers[i].barrier1_time=0;
    workers[i].barrier2_time=0;
    workers[i].verticesProcessed=0;
    workers[i].edgesProcessed=0;
    workers[i].nextVertexToBeProcessed_time=0;
  }
  for(int i=0;i<g.n_;i++)
  {
    atomics_current[i].store(pr_curr[i]);
    atomics_next[i].store(pr_next[i]);
  }
  timer t1;
  double time_taken=0.0;
  t1.start();
  for(int i=0;i<n_threads;i++)
  {
    threads[i]= std::thread(Distributed_PageRank_Calculation,&workers[i]);
  }
  for(int i=0;i<n_threads;i++)
  {
    threads[i].join();
  }
  time_taken= t1.stop();
  for(int i=0;i<g.n_;i++)
  {
    pr_curr[i]=atomics_current[i].load();
  }


  PageRankType sum_of_page_ranks = 0;

  for (uintV u = 0; u < g.n_; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
 
  
 
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time" << std::endl;
  for(int i=0;i<n_threads;i++)
  {
    std::cout<<i<<", "<<workers[i].verticesProcessed<<", "<<workers[i].edgesProcessed
    <<", "<<barrier1_time[i]<<std::setprecision(6)<<", "<<barrier2_time[i]<<std::setprecision(6)<<", "<<
    nextVertexToBeProcessed_time[i]<<std::setprecision(6)<<", "<<timeTaken[i]<<std::endl;
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks <<std::setprecision(6)<<"\n";
  std::cout << "Time taken (in seconds) : " << time_taken << std::setprecision(6)<<"\n";
  
  pageRank_Destroy();
  delete[] timeTaken;
  delete[] barrier1_time;
  delete[] barrier2_time;
  delete[] nextVertexToBeProcessed_time; 
  delete[] atomics_current;
  delete[] atomics_next;
  delete barrier;
  
}
void strategy2(uint &n_threads, Graph&g, uint &max_iterations)
{
  pageRank_Init(g);
  CustomBarrier* barrier= new CustomBarrier(n_threads);    
  timeTaken= new double[n_threads]();
  barrier1_time= new double[n_threads]();
  barrier2_time= new double[n_threads]();
  nextVertexToBeProcessed_time= new double[n_threads]();
  // Initializing worker threads
  std::thread threads[n_threads];
  worker workers[n_threads];
  
  /*
  *Distributing vertices between threads based on its outdegree
  */  
  std::vector<uint> startx(n_threads);
  std::vector<uint> endx(n_threads);
  
  uint numOfEdgesForThread= g.m_/n_threads;
  uint totalAssignedEdges=0;
  int currentlyProcessingVertex=0;
  for(int i=0;i< n_threads;i++){

    startx[i]=currentlyProcessingVertex;
    while (totalAssignedEdges<(i+1)*numOfEdgesForThread)
    {
      //get out degree of a vertex
      uint outDegree= g.vertices_[currentlyProcessingVertex].getOutDegree();
      currentlyProcessingVertex++;
      totalAssignedEdges+=outDegree;
      endx[i]=currentlyProcessingVertex-1;
    }
    
    if(currentlyProcessingVertex==g.n_-1){
      endx[i]=currentlyProcessingVertex;
    }

  } 

  atomics_current= new std::atomic<PageRankType>[g.n_];
  atomics_next= new std::atomic<PageRankType>[g.n_];

  for(int i=0;i< n_threads;i++)
  {
    workers[i].barrier= barrier;
    workers[i].startVertex= startx[i];
    workers[i].endVertex= endx[i];
    workers[i].graph_ptr= &g;
    workers[i].pr_curr= pr_curr;
    workers[i].pr_next= pr_next;
    workers[i].threadID=i;
    workers[i].time_buffer=timeTaken;
    workers[i].total_num_vertices=g.n_;
    workers[i].max_iters= max_iterations;
    workers[i].atomics= atomics_current;
    workers[i].barrier1_time=0;
    workers[i].barrier2_time=0;
    workers[i].verticesProcessed=0;
    workers[i].edgesProcessed=0;
    workers[i].nextVertexToBeProcessed_time=0;
 
  }
 
  for(int i=0;i<g.n_;i++)
  {
    atomics_current[i].store(pr_curr[i]);
    atomics_next[i].store(pr_next[i]);
  }
  timer t1;
  double time_taken=0.0;
  t1.start();
  for(int i=0;i<n_threads;i++)
  {
    threads[i]= std::thread(Distributed_PageRank_Calculation,&workers[i]);
  }
  for(int i=0;i<n_threads;i++)
  {
    threads[i].join();
  }
  time_taken= t1.stop();
  for(int i=0;i<g.n_;i++)
  {
    pr_curr[i]=atomics_current[i].load();
  }


  PageRankType sum_of_page_ranks = 0;

  for (uintV u = 0; u < g.n_; u++) {
    sum_of_page_ranks += pr_curr[u];
  }

  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time" << std::endl;
  for(int i=0;i<n_threads;i++)
  {
    std::cout<<i<<", "<<workers[i].verticesProcessed<<", "<<workers[i].edgesProcessed
    <<", "<<barrier1_time[i]<<std::setprecision(6)<<", "<<barrier2_time[i]<<std::setprecision(6)<<", "<<
    nextVertexToBeProcessed_time[i]<<std::setprecision(6)<<", "<<timeTaken[i]<<std::endl;
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks <<std::setprecision(6)<<"\n";
  std::cout << "Time taken (in seconds) : " << time_taken << std::setprecision(6)<<"\n";
  /*
  * Clean-up 
  */
  pageRank_Destroy();
  delete[] timeTaken;
  delete[] barrier1_time;
  delete[] barrier2_time;
  delete[] nextVertexToBeProcessed_time; 
  delete[] atomics_current;
  delete[] atomics_next;
  delete barrier;
}
static std::atomic<int> nextVertex(0);
//static std::atomic<uint> numberOfCalls(0);
int getNextVertexToBeProcessed(uint graph_size,uint granularity)
{
 
  int old= nextVertex.load(); 
  int next=old+granularity;
  while(!nextVertex.compare_exchange_weak(old, next))
  {
     next=old+granularity;
  }
  if(old>=graph_size)
  {
    return -1;
  }

  //numberOfCalls++;
  return old;
}
inline void resetNextVertexToBeProcessed()
{
  nextVertex=0;
}
void DynamicMapping_PageRankPull_Calculation(worker* worker)
{
  timer t;
  timer barrier1_timer;
  timer barrier2_timer; 
  timer nextVertex_timer; 
  double time_taken=0.0;
  double barrier1_time_taken=0.0;
  double barrier2_time_taken=0.0;
  double nextVertex_time_taken=0.0; 
  PageRankType* prr_curr_local= worker->pr_curr;
  PageRankType* prr_next_local= worker->pr_next;  
  uint max_iters= worker->max_iters;
  CustomBarrier *barrier_local= worker->barrier;
  Graph* g= worker->graph_ptr;
  std:: atomic <PageRankType>* atomics_local= worker->atomics;

  uint edgesProcessed=0;
  uint verticesProcessed=0;
  t.start();
  for (int iter = 0; iter < max_iters; iter++) {
    
    //---------------------------------------------------------
    while(1)
    {
      nextVertex_timer.start();
      int u= getNextVertexToBeProcessed(g->n_,1);
      nextVertex_time_taken+=nextVertex_timer.stop();
      if(u==-1)
      {
        break;
      }  
      uintE out_degree = g->vertices_[u].getOutDegree();
      edgesProcessed+=out_degree;
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g->vertices_[u].getOutNeighbor(i);
        
        //compare and swap here
        PageRankType fetch_old = atomics_next[v].load();

        PageRankType sum = fetch_old+(atomics_current[u].load() / (PageRankType) out_degree);

        while(!atomics_next[v].compare_exchange_weak(fetch_old, sum))
        {
          sum= fetch_old+(atomics_current[u].load() / (PageRankType) out_degree);
        }
        
      }
    }
    barrier1_timer.start();
    barrier_local->wait();
    barrier1_time_taken+= barrier1_timer.stop();
    /*
    Reset the next vertex to be processed to 0 if it's thread 0
    */
    //-----------------------------------------------------------
    if(worker->threadID==0)
    {
      resetNextVertexToBeProcessed();
      barrier_local->wait();
    }
    else
    {
      barrier_local->wait();
    }
    //-----------------------------------------------------------
    while(1)
    {
      nextVertex_timer.start();
      int v= getNextVertexToBeProcessed(g->n_,1);   
      nextVertex_time_taken+=nextVertex_timer.stop();
      if(v==-1){
        break;
      } 
      verticesProcessed++;
      atomics_next[v].store(PAGE_RANK(atomics_next[v].load()));
      // reset pr_curr for the next iteration
      atomics_current[v].store(atomics_next[v].load());
      atomics_next[v].store (0.0);
    }
    barrier2_timer.start();
    barrier_local->wait();
    barrier2_time_taken+=barrier2_timer.stop();
    /*
    Reset the next vertex to be processed to 0 if it's thread 0
    */
    //----------------------------------------------------------
    if(worker->threadID==0)
    {
      resetNextVertexToBeProcessed();
      barrier_local->wait();
    }
    else
    {
      barrier_local->wait();
    }   
  }
  time_taken=t.stop();
  worker->time_buffer[worker->threadID]+=time_taken;
  worker->edgesProcessed=edgesProcessed;
  worker->verticesProcessed=verticesProcessed;
  worker->nextVertexToBeProcessed_time[worker->threadID]+= nextVertex_time_taken;
  worker->barrier1_time[worker->threadID]+=barrier1_time_taken;
  worker->barrier2_time[worker->threadID]+=barrier2_time_taken;
}
void strategy3(uint &n_threads, Graph&g, uint &max_iterations)
{
  pageRank_Init(g);
  CustomBarrier* barrier= new CustomBarrier(n_threads);    
  timeTaken= new double[n_threads]();
  barrier1_time= new double[n_threads]();
  barrier2_time= new double[n_threads]();
  nextVertexToBeProcessed_time= new double[n_threads]();
  // Initializing worker threads
  std::thread threads[n_threads];
  worker workers[n_threads];

  for(int i=0;i< n_threads;i++)
  {
    workers[i].barrier= barrier;
    workers[i].graph_ptr= &g;
    workers[i].pr_curr= pr_curr;
    workers[i].pr_next= pr_next;
    workers[i].threadID=i;
    workers[i].time_buffer=timeTaken;
    workers[i].total_num_vertices=g.n_;
    workers[i].max_iters= max_iterations;
    workers[i].atomics= atomics_current;
    workers[i].barrier1_time=barrier1_time;
    workers[i].barrier2_time=barrier2_time;
    workers[i].nextVertexToBeProcessed_time= nextVertexToBeProcessed_time;
   
  }
  atomics_current= new std::atomic<PageRankType>[g.n_];
  atomics_next= new std::atomic<PageRankType>[g.n_];
  for(int i=0;i<g.n_;i++)
  {
    atomics_current[i].store(pr_curr[i]);
    atomics_next[i].store(pr_next[i]);
  }
  timer t1;
  double time_taken=0.0;
  t1.start();
  for(int i=0;i<n_threads;i++)
  {
    threads[i]= std::thread(DynamicMapping_PageRankPull_Calculation,&workers[i]);
  }
  for(int i=0;i<n_threads;i++)
  {
    threads[i].join();
  }
  time_taken= t1.stop();
  for(int i=0;i<g.n_;i++)
  {
    pr_curr[i]=atomics_current[i].load();
  }


  PageRankType sum_of_page_ranks = 0;

  for (uintV u = 0; u < g.n_; u++) {
    sum_of_page_ranks += pr_curr[u];
  }

  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time" << std::endl;
  for(int i=0;i<n_threads;i++)
  {
    std::cout<<i<<", "<<workers[i].verticesProcessed<<", "<<workers[i].edgesProcessed
    <<", "<<barrier1_time[i]<<std::setprecision(6)<<", "<<barrier2_time[i]<<std::setprecision(6)<<", "<<
    nextVertexToBeProcessed_time[i]<<std::setprecision(6)<<", "<<timeTaken[i]<<std::endl;
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks <<std::setprecision(6)<<"\n";
  std::cout << "Time taken (in seconds) : " << time_taken << std::setprecision(6)<<"\n";
  
  /*
  * Clean-up 
  */
  pageRank_Destroy();
  delete[] timeTaken;
  delete[] barrier1_time;
  delete[] barrier2_time;
  delete[] nextVertexToBeProcessed_time; 
  delete[] atomics_current;
  delete[] atomics_next;
  delete barrier;
}

void CoarseGrained_DynamicMapping_PageRankPull_Calculation(worker* worker)
{
  timer t;
  timer barrier1_timer;
  timer barrier2_timer; 
  timer nextVertex_timer; 
  double time_taken=0.0;
  double barrier1_time_taken=0.0;
  double barrier2_time_taken=0.0;
  double nextVertex_time_taken=0.0; 
  PageRankType* prr_curr_local= worker->pr_curr;
  PageRankType* prr_next_local= worker->pr_next;  
  uint max_iters= worker->max_iters;
  CustomBarrier *barrier_local= worker->barrier;
  Graph* g= worker->graph_ptr;
  std:: atomic <PageRankType>* atomics_local= worker->atomics;
  uint k= worker->granularity;
  uint edgesProcessed=0;
  uint verticesProcessed=0;
  t.start();
  for (int iter = 0; iter < max_iters; iter++) {
    
    while(1)
    {
      nextVertex_timer.start();
      int u= getNextVertexToBeProcessed(g->n_, k);
      nextVertex_time_taken+=nextVertex_timer.stop();

     // std::cout<<"U= "<<u<<" Thread# "<<worker->threadID<<std::endl;
      if(u==-1)
      {
        break;
      }  
      
      for(int j=0;j<k;j++)
      {
        uintE out_degree = g->vertices_[u].getOutDegree();
        edgesProcessed+=out_degree;
        for (uintE i = 0; i < out_degree; i++) {
          uintV v = g->vertices_[u].getOutNeighbor(i);
          
          //compare and swap here
          PageRankType fetch_old = atomics_next[v].load();

          PageRankType sum = fetch_old+(atomics_current[u].load() / (PageRankType) out_degree);

          while(!atomics_next[v].compare_exchange_weak(fetch_old, sum))
          {
            sum= fetch_old+(atomics_current[u].load() / (PageRankType) out_degree);
          }
          
        }
        u++;
        if(u>=g->n_)
        {
          break;
        }
      }
    }
    barrier1_timer.start();
    barrier_local->wait();
    barrier1_time_taken+= barrier1_timer.stop();
  
    /*
    Reset the number of processed vertices to 0 if it's thread 0
    */
    //-----------------------------------------------------------
    if(worker->threadID==0)
    {
      resetNextVertexToBeProcessed();
      barrier_local->wait();
    }
    else
    {
      barrier_local->wait();
    }
    //-----------------------------------------------------------
    
    while(1)
    {
      nextVertex_timer.start();
      int v= getNextVertexToBeProcessed(g->n_,k);      
      nextVertex_time_taken+=nextVertex_timer.stop();
      if(v==-1){
        break;
      } 
      for(int j=0;j<k;j++){
        verticesProcessed++;
        atomics_next[v].store(PAGE_RANK(atomics_next[v].load()));
        // reset pr_curr for the next iteration
        atomics_current[v].store(atomics_next[v].load());
        atomics_next[v].store (0.0);
        v++;
        if(v>=g->n_)
        {
          break;
        }
      }
    }
    barrier2_timer.start();
    barrier_local->wait();
    barrier2_time_taken+=barrier2_timer.stop();
    /*
    Reset the number of processed vertices to 0 if it's thread 0
    */
    //----------------------------------------------------------
    if(worker->threadID==0)
    {
      resetNextVertexToBeProcessed();
      barrier_local->wait();
    }
    else
    {
      barrier_local->wait();
    }
    //---------------------------------------------------------
       
  }
  time_taken=t.stop();
  worker->time_buffer[worker->threadID]+=time_taken;
  worker->edgesProcessed=edgesProcessed;
  worker->verticesProcessed=verticesProcessed;
  worker->nextVertexToBeProcessed_time[worker->threadID]+= nextVertex_time_taken;
  worker->barrier1_time[worker->threadID]+=barrier1_time_taken;
  worker->barrier2_time[worker->threadID]+=barrier2_time_taken;
}
void strategy4(uint &n_threads, Graph&g, uint &max_iterations, uint &granularity)
{
  pageRank_Init(g);
  CustomBarrier* barrier= new CustomBarrier(n_threads);    
  timeTaken= new double[n_threads]();
  barrier1_time= new double[n_threads]();
  barrier2_time= new double[n_threads]();
  nextVertexToBeProcessed_time= new double[n_threads]();
  
  // Initializing worker threads
  std::thread threads[n_threads];
  worker workers[n_threads];

  for(int i=0;i< n_threads;i++)
  {
    workers[i].barrier= barrier;
    workers[i].graph_ptr= &g;
    workers[i].pr_curr= pr_curr;
    workers[i].pr_next= pr_next;
    workers[i].threadID=i;
    workers[i].time_buffer=timeTaken;
    workers[i].total_num_vertices=g.n_;
    workers[i].max_iters= max_iterations;
    workers[i].atomics= atomics_current;
    workers[i].granularity=granularity;
    workers[i].barrier1_time=barrier1_time;
    workers[i].barrier2_time=barrier2_time;
    workers[i].nextVertexToBeProcessed_time= nextVertexToBeProcessed_time;
  }
  atomics_current= new std::atomic<PageRankType>[g.n_];
  atomics_next= new std::atomic<PageRankType>[g.n_];
  for(int i=0;i<g.n_;i++)
  {
    atomics_current[i].store(pr_curr[i]);
    atomics_next[i].store(pr_next[i]);
  }
  timer t1;
  double time_taken=0.0;
  t1.start();
  for(int i=0;i<n_threads;i++)
  {
    threads[i]= std::thread(CoarseGrained_DynamicMapping_PageRankPull_Calculation,&workers[i]);
  }
  for(int i=0;i<n_threads;i++)
  {
    threads[i].join();
  }
  time_taken= t1.stop();
  for(int i=0;i<g.n_;i++)
  {
    pr_curr[i]=atomics_current[i].load();
  }


  PageRankType sum_of_page_ranks = 0;

  for (uintV u = 0; u < g.n_; u++) {
    sum_of_page_ranks += pr_curr[u];
  }

  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time" << std::endl;
  for(int i=0;i<n_threads;i++)
  {
    std::cout<<i<<", "<<workers[i].verticesProcessed<<", "<<workers[i].edgesProcessed
    <<", "<<barrier1_time[i]<<std::setprecision(6)<<", "<<barrier2_time[i]<<std::setprecision(6)<<", "<<
    nextVertexToBeProcessed_time[i]<<std::setprecision(6)<<", "<<timeTaken[i]<<std::endl;
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks <<std::setprecision(6)<<"\n";
  std::cout << "Time taken (in seconds) : " << time_taken << std::setprecision(6)<<"\n";
  
  /*
  * Clean-up 
  */
  pageRank_Destroy();
  delete[] timeTaken;
  delete[] barrier1_time;
  delete[] barrier2_time;
  delete[] nextVertexToBeProcessed_time; 
  delete[] atomics_current;
  delete[] atomics_next;
  delete barrier;
}
int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nThreads", "Number of Threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
          {"strategy", "Task decomposition and mapping strategy number",
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
          {"granularity", "Granulity for strategy4",
           cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  uint strategy= cl_options["strategy"].as<uint>();
  uint granularity= cl_options["granularity"].as<uint>();
#ifdef USE_INT
  std::cout << "Using INT" << std::endl;
#else
  std::cout << "Using DOUBLE" << std::endl;
#endif
  if(strategy!=4)
  {
    granularity=1;
  }
  if(granularity<1 && strategy==4)
  {
    std::cout<<"Invalid granularity argument, terminating..."<<std::endl;
    exit(-1);
  }
  std::cout << std::fixed;
  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout<<"Strategy : "<<strategy<<std::endl;
  std::cout<<"Granularity : "<<granularity<<std::endl;
  std::cout << "Iterations: " <<max_iterations << std::endl;
  
  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  switch (strategy)
  {
  case 1:
    strategy1(n_threads, g, max_iterations);
    break;
  
  case 2:
    strategy2(n_threads,g, max_iterations);
    break;
  case 3:
    strategy3(n_threads, g,max_iterations);
    break;
  case 4:
    strategy4(n_threads,g, max_iterations, granularity);
    break;
  default:
    std::cout<<"Unknown strategy"<<std::endl;
    break;
  }
  return 0;
}
