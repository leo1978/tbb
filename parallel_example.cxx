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

class GenerateContour
{
  vtkImageData** Images;

public:
  void operator() ( const blocked_range<size_t>& r ) const
    {
      vtkContourFilter* cf = vtkContourFilter::New();
      cf->SetValue(0, 200);
      vtkImageData** images = this->Images;
      for ( size_t i = r.begin(); i != r.end(); ++i )
        {
        vtkImageData* image = images[i];
        cf->SetInputData(image);
        cf->Update();
        if (!cf->GetOutput()->GetNumberOfCells())
          abort();
        }
      cf->Delete();
    }
  GenerateContour(vtkImageData** s) : Images(s)
    {
    }
};

class GenerateContour2
{
  vtkImageData** Images;

public:
  void operator() ( const blocked_range<size_t>& r ) const
    {
      vtkImageData** images = this->Images;

      vtkContourFilter* cf = vtkContourFilter::New();
      cf->SetValue(0, 200);
      cf->SetInputData(images[r.begin()]);
      cf->PropagateUpdateExtent();

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
      cf->Delete();
    }
  GenerateContour2(vtkImageData** s) : Images(s)
    {
    }
};

class GenerateContour3
{
  vtkImageData** Images;

public:
  void operator() ( const blocked_range<size_t>& r ) const
    {
      vtkImageData** images = this->Images;

      vtkContourFilter* cf = vtkContourFilter::New();
      cf->SetValue(0, 200);
      cf->SetInputData(images[r.begin()]);
      cf->PropagateUpdateExtent(); // Create output

      vtkInformation* request1 = vtkInformation::New();
      request1->Set(vtkDemandDrivenPipeline::REQUEST_INFORMATION());

      vtkInformation* request2 = vtkInformation::New();
      request2->Set(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT());

      vtkInformation* request3 = vtkInformation::New();
      request3->Set(vtkDemandDrivenPipeline::REQUEST_DATA());

      vtkInformation* inInfo = cf->GetInputInformation();
      vtkInformation* outInfo = cf->GetOutputInformation(0);
      vtkExecutive* executive = cf->GetExecutive();

      for ( size_t i = r.begin(); i != r.end(); ++i )
        {
        vtkImageData* image = images[i];
        inInfo->Set(vtkDataObject::DATA_OBJECT(), image);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
                    image->GetExtent(),
                    6);
        cf->ProcessRequest(request1,
                           executive->GetInputInformation(),
                           executive->GetOutputInformation());
        outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                     image->GetExtent(),
                     6);
        cf->ProcessRequest(request2,
                           executive->GetInputInformation(),
                           executive->GetOutputInformation());
        cf->ProcessRequest(request3,
                           executive->GetInputInformation(),
                           executive->GetOutputInformation());
        if (!cf->GetOutput()->GetNumberOfCells())
          abort();
        }
      cf->Delete();
    }
  GenerateContour3(vtkImageData** s) : Images(s)
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

  GenerateContour2 gc2( images );

  GenerateContour3 gc3( images );

  vtkTimerLog::MarkStartEvent("serial");
  gc(range);
  vtkTimerLog::MarkEndEvent("serial");

  vtkTimerLog::MarkStartEvent("serial2");
  gc2(range);
  vtkTimerLog::MarkEndEvent("serial2");

  vtkTimerLog::MarkStartEvent("serial3");
  gc3(range);
  vtkTimerLog::MarkEndEvent("serial3");

  vtkTimerLog::MarkStartEvent("parallel");
  parallel_for(blocked_range<size_t>(0, N, 20), gc);
  vtkTimerLog::MarkEndEvent("parallel");

  vtkTimerLog::MarkStartEvent("parallel2");
  parallel_for(blocked_range<size_t>(0, N, 20), gc2);
  vtkTimerLog::MarkEndEvent("parallel2");

  vtkTimerLog::MarkStartEvent("parallel3");
  parallel_for(blocked_range<size_t>(0, N, 20), gc3);
  vtkTimerLog::MarkEndEvent("parallel3");

  vtkTimerLog::DumpLogWithIndents(&std::cout, 0.00001);

  source->Delete();
  for(int i=0; i<N; i++)
    {
    images[i]->Delete();
    }

  return 0;
}

