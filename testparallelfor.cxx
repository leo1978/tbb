#include "vtkParallelFor.h"
#include "vtkTimerLog.h"

#include "tbb/tbb.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include <tbb/atomic.h>
#include "tbb/mutex.h"
#include "assert.h"

#include <pthread.h>
#include <iostream>
using namespace std;
using namespace tbb;

int TotalWork = 0;
int MinWork=INT_MAX;
int MaxWork=0;
bool SingleThread = true;
class LoopLots
{
public:
  LoopLots(int m, int n, double k): M(m), N(n), K(k)
  {
    assert(m<=n);
  }

  void operator()(const blocked_range<size_t>& r) const
  {
    double k = 1.0;
    // PrintLock.lock();
    // cout<<r.begin()<<" "<<r.end()<<endl;
    // PrintLock.unlock();
    for(size_t i = r.begin(); i!=r.end(); i++)
      {
      int work = static_cast<int>(K*(i - M/2.0) + N/2.0);
      assert(work>=0);
      if(SingleThread)
        {
        TotalWork+=work;
        MinWork = min(MinWork,work);
        MaxWork = max(MaxWork,work);
        }
      int s =0;
      for(int ii=0; ii<work; ii++)
        {
        for(int jj=0; jj<100; jj++)
          {
          s = s+1;
          }
        }
      }
  }
private:
  int M;
  int N;
  double K;
};

int main(int argc, char* argv[])
{
  if(argc==1)
    {
    cout<<"Run with options: [balance] [-vtk|-tbb|-1]"<<endl;
    return 0;
    }
  bool useVTK = argc==3 && strcmp(argv[1],"-vtk");
  int work(10000);
  blocked_range<size_t> range(0,work);
  SingleThread = argc==3 && strcmp(argv[2],"-1")==0;
  double unbalance =  atof(argv[1]);
  LoopLots f(work,work,unbalance);

  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  if(SingleThread)
    {
    cout<<"Single Thread Execute\n";
    f(range);
    }
  else if (useVTK)
    {
    cout<<"vtkParallelFor execute\n";
    vtkParallelFor(range,f);
    }
  else
    {
    cout<<"tbb execute\n";
    parallel_for(range, f);
    }
  timer->StopTimer();
  std::cout << "time: " << timer->GetElapsedTime() <<" seconds"<< std::endl;

  if(SingleThread)
    {
    cout<<"stats: "<<TotalWork<<" "<<MinWork<<" "<<MaxWork<<endl;
    }
  return 0;
}
