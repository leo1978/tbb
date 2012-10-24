#include <iostream>
#include <string>
#include <algorithm>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <vtkSynchronizedTemplates3D.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkRTAnalyticSource.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTimerLog.h>

using namespace tbb;
using namespace std;

static const size_t N = 2000;
static const size_t GRAIN = 20;

class GenerateContour
{
  vtkImageData** Images;
  vtkPolyData** Outputs;
  vtkSynchronizedTemplates3D* ContourFilter;

public:
  void operator() ( const blocked_range<size_t>& r ) const
    {
      vtkImageData** images = this->Images;
      vtkPolyData** outputs = this->Outputs;
      vtkSynchronizedTemplates3D* cf = this->ContourFilter;

      vtkInformation* request1 = vtkInformation::New();
      request1->Set(vtkDemandDrivenPipeline::REQUEST_INFORMATION());

      vtkInformation* request2 = vtkInformation::New();
      request2->Set(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT());

      vtkInformation* request3 = vtkInformation::New();
      request3->Set(vtkDemandDrivenPipeline::REQUEST_DATA());

      vtkInformationVector* inInfoVec = vtkInformationVector::New();
      vtkInformation* inInfo = vtkInformation::New();
      inInfoVec->Append(inInfo);
      inInfo->Delete();

      vtkInformationVector* outInfoVec = vtkInformationVector::New();
      vtkInformation* outInfo = vtkInformation::New();
      outInfoVec->Append(outInfo);
      outInfo->Delete();

      vtkExecutive* executive = cf->GetExecutive();

      for ( size_t i = r.begin(); i != r.end(); ++i )
        {
        vtkImageData* image = images[i];
        vtkPolyData* output = outputs[i];
        //vtkPolyData* output = vtkPolyData::New();
        outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
        inInfo->Set(vtkDataObject::DATA_OBJECT(), image);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
                    image->GetExtent(),
                    6);
        cf->ProcessRequest(request1,
                           &inInfoVec,
                           outInfoVec);

        outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                     image->GetExtent(),
                     6);
        cf->ProcessRequest(request2,
                           &inInfoVec,
                           outInfoVec);

        cf->ProcessRequest(request3,
                           &inInfoVec,
                           outInfoVec);
        //output->Initialize();
        /*
        if (!output->GetNumberOfCells())
          abort();
        */
        //outputs[i]->ShallowCopy(output);
        //output->Delete();
        }
      inInfoVec->Delete();
      outInfoVec->Delete();
    }
  GenerateContour(vtkImageData** s, vtkPolyData** opts,
                  vtkSynchronizedTemplates3D* cf) :
    Images(s), Outputs(opts), ContourFilter(cf)
    {
    }
};

int main()
{
  vtkRTAnalyticSource* source = vtkRTAnalyticSource::New();
  //source->SetWholeExtent(0, 40, 0, 40, 0, 40);
  //source->SetWholeExtent(0, 80, 0, 80, 0, 80);
  source->Update();
  std::cout << "Number of cells per block: "
            << source->GetOutput()->GetNumberOfCells()
            << std::endl;

  vtkImageData* images[N];
  vtkPolyData* pd[N];
  for(int i=0; i<N; i++)
    {
    images[i] = vtkImageData::New();
    images[i]->DeepCopy(source->GetOutput());
    pd[i] = vtkPolyData::New();
    }

  vtkSynchronizedTemplates3D* cf = vtkSynchronizedTemplates3D::New();
  cf->SetValue(0, 200);

  GenerateContour gc( images, pd, cf );
  blocked_range<size_t> range(0, N);

  vtkTimerLog::MarkStartEvent("serial");
  gc(range);
  vtkTimerLog::MarkEndEvent("serial");

  vtkTimerLog::MarkStartEvent("parallel");
  parallel_for(blocked_range<size_t>(0, N, GRAIN), gc);
  vtkTimerLog::MarkEndEvent("parallel");

  vtkTimerLog::DumpLogWithIndents(&std::cout, 0.00001);

  cf->Delete();
  source->Delete();
  for(int i=0; i<N; i++)
    {
    images[i]->Delete();
    pd[i]->Delete();
    }

  return 0;
}

