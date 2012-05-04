#include <iostream>
#include <string>
#include <algorithm>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkRTAnalyticSource.h>
#include <vtkTimerLog.h>

using namespace tbb;
using namespace std;
static const size_t N = 200;

class GenerateContour
{
  vtkImageData** Images;

public:
  void operator() ( const blocked_range<size_t>& r ) const 
    {
      vtkContourFilter* cf = vtkContourFilter::New();
      vtkImageData** images = this->Images;
      for ( size_t i = r.begin(); i != r.end(); ++i )
        {
        vtkImageData* image = images[i];
        cf->SetInputData(image);
        cf->Update();
        }
      cf->Delete();
    }
  GenerateContour(vtkImageData** s) : Images(s) 
    {
    }
};

int main()
{
  vtkRTAnalyticSource* source = vtkRTAnalyticSource::New();
  //source->SetWholeExtent(0, 80, 0, 80, 0, 80);
  source->Update();
  std::cout << "Number of cells per block: "
            << source->GetOutput()->GetNumberOfCells()
            << std::endl;

  vtkImageData* images[N];
  for(int i=0; i<N; i++)
    {
    images[i] = vtkImageData::New();
    images[i]->DeepCopy(source->GetOutput());
    }

  GenerateContour gc( images );
  blocked_range<size_t> range(0, N);

  vtkTimerLog::MarkStartEvent("serial");
  gc(range);
  vtkTimerLog::MarkEndEvent("serial");

  vtkTimerLog::MarkStartEvent("parallel");
  parallel_for(blocked_range<size_t>(0, N, 20), gc);
  vtkTimerLog::MarkEndEvent("parallel");

  vtkTimerLog::DumpLogWithIndents(&std::cout, 0.00001);

  source->Delete();
  for(int i=0; i<N; i++)
    {
    images[i]->Delete();
    }

  return 0;
}

