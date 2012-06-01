#include <iostream>
#include <string>
#include <algorithm>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <vtkContourFilter.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkRTAnalyticSource.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTimerLog.h>

using namespace tbb;
using namespace std;

static const size_t N = 2000;
static const size_t GRAIN = 1;

class TestThreadSafety
{
  vtkImageData** Images;

public:
  void operator() ( const blocked_range<size_t>& r ) const
    {
      vtkImageData** images = this->Images;

      vtkContourFilter* cf = vtkContourFilter::New();
      cf->SetValue(0, 200);
      cf->SetInputData(images[0]);
      cf->PropagateUpdateExtent();

      /*
      vtkInformation* request = vtkInformation::New();
      request->Set(vtkDemandDrivenPipeline::REQUEST_DATA());

      vtkInformation* inInfo = cf->GetInputInformation();
      vtkExecutive* executive = cf->GetExecutive();

      for ( size_t i = r.begin(); i != r.end(); ++i )
        {
        vtkImageData* image = images[i];
        inInfo->Set(vtkDataObject::DATA_OBJECT(), image);
        cf->ProcessRequest(request,
                           executive->GetInputInformation(),
                           executive->GetOutputInformation());

        if (!cf->GetOutput()->GetNumberOfCells())
          abort();
        }
      */

      cf->Delete();
    }
  TestThreadSafety(vtkImageData** s) : Images(s)
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

  TestThreadSafety test( images );

  vtkTimerLog::MarkStartEvent("parallel");
  parallel_for(blocked_range<size_t>(0, N, GRAIN), test);
  vtkTimerLog::MarkEndEvent("parallel");

  std::cout << images[0]->GetReferenceCount() << endl;

  source->Delete();
  for(int i=0; i<N; i++)
    {
    images[i]->Delete();
    if (i == 0)
      std::cout << images[0]->GetReferenceCount() << endl;
    }


  return 0;
}

