#include <iostream>

#include "tbb/atomic.h"

#include <vtkCriticalSection.h>
#include <vtkMultiThreader.h>

using namespace tbb;
using namespace std;

// This program demonstrates the use of the atomic operation fetch_and_add().
// It also empirically demonstrates that it works. I use the following to test
// this:
// for i in `seq 10000`
// do
// ./atomic_example | sort -g | uniq | wc -l
// done | grep -v 14
//
// Run the script above for both cases. If thread safe, it should not produce
// any output.
atomic<unsigned int> counter;

const unsigned int NUM_THREADS = 14;
unsigned int vals[NUM_THREADS];

double foo = 0;

void* Threaded(void* arg)
{
  vtkMultiThreader::ThreadInfo* tInfo =
    (vtkMultiThreader::ThreadInfo*)arg;

  // Not thread safe
  //counter += 1;
  //vals[tInfo->ThreadID] = counter;

  // Thread safe
  vals[tInfo->ThreadID] = counter.fetch_and_add(1) + 1;

  return 0;
}

int main()
{
  counter = 0;

  vtkMultiThreader* mt = vtkMultiThreader::New();
  mt->SetSingleMethod(Threaded, 0);
  mt->SetNumberOfThreads(NUM_THREADS);
  mt->SingleMethodExecute();

  for(unsigned int i=0; i<NUM_THREADS; i++)
    cout << vals[i] << endl;

  mt->Delete();
  return 0;
}

