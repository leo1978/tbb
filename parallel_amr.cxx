#include <iostream>
#include <string>
#include <algorithm>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <vtkSynchronizedTemplates3D.h>
#include <vtkAMRFlashReader.h>
#include <vtkContourFilter.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTimerLog.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>

#include <vtkOverlappingAMR.h>
#include <vtkCompositeDataIterator.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>

using namespace tbb;
using namespace std;

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

      vtkNew<vtkCellDataToPointData> c2p;

      for ( size_t i = r.begin(); i != r.end(); ++i )
        {
        vtkImageData* image = images[i];

        c2p->SetInputData(image);
        c2p->Update();

        image = dynamic_cast<vtkImageData*>(c2p->GetOutput());

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
  vtkNew<vtkAMRFlashReader> reader;
  reader->SetFileName("/Users/berk/Desktop/smooth/smooth.flash");
  reader->SetMaxLevel(7);
  reader->SetCellArrayStatus("dens", 1);
  reader->Update();

  std::vector<vtkImageData*> images;
  vtkSmartPointer<vtkCompositeDataIterator> iter;
  iter.TakeReference(reader->GetOutput()->NewIterator());

  cout << "here" << endl;

  for(iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
    vtkImageData* id = dynamic_cast<vtkImageData*>(iter->GetCurrentDataObject());
    if (id)
      {
      images.push_back(id);
      double* range = id->GetCellData()->GetArray("dens")->GetRange(0);
      }
    }

  std::vector<vtkPolyData*> outputs(images.size());
  std::vector<vtkPolyData*>::iterator iter2;
  for (iter2 = outputs.begin();
       iter2 != outputs.end();
       iter2++)
    {
    *iter2 = vtkPolyData::New();
    }

  vtkNew<vtkSynchronizedTemplates3D> cf;
  cf->SetInputArrayToProcess(0, 0, 0, 0, "dens");
  cf->SetValue(0, 3000);

  GenerateContour gc( &images[0], &outputs[0], cf.GetPointer() );
  blocked_range<size_t> range(0, images.size());

  vtkTimerLog::MarkStartEvent("serial");
  gc(range);
  vtkTimerLog::MarkEndEvent("serial");

  for (iter2 = outputs.begin();
       iter2 != outputs.end();
       iter2++)
    {
    (*iter2)->Delete();
    *iter2 = vtkPolyData::New();
    }

  vtkTimerLog::MarkStartEvent("parallel");
  parallel_for(blocked_range<size_t>(0, images.size(), GRAIN), gc);
  vtkTimerLog::MarkEndEvent("parallel");

  vtkTimerLog::DumpLogWithIndents(&std::cout, 0.00001);

  vtkNew<vtkMultiBlockDataSet> mb;
  mb->SetNumberOfBlocks(outputs.size());
  int i;
  for (iter2 = outputs.begin(), i=0;
       iter2 != outputs.end();
       iter2++, i++)
    {
    mb->SetBlock(i, *iter2);
    }

  vtkNew<vtkXMLMultiBlockDataWriter> writer;
  writer->SetInputData(mb.GetPointer());
  writer->SetFileName("/Users/berk/polys.vtm");
  //writer->Write();

  for (iter2 = outputs.begin();
       iter2 != outputs.end();
       iter2++)
    {
    (*iter2)->Delete();
    }

  /*
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
  parallel_for(blocked_range<size_t>(0, N, GRAIN), gc);
  vtkTimerLog::MarkEndEvent("parallel");

  vtkTimerLog::MarkStartEvent("parallel2");
  parallel_for(blocked_range<size_t>(0, N, GRAIN), gc2);
  vtkTimerLog::MarkEndEvent("parallel2");

  vtkTimerLog::MarkStartEvent("parallel3");
  parallel_for(blocked_range<size_t>(0, N, GRAIN), gc3);
  vtkTimerLog::MarkEndEvent("parallel3");

  vtkTimerLog::DumpLogWithIndents(&std::cout, 0.00001);

  source->Delete();
  for(int i=0; i<N; i++)
    {
    images[i]->Delete();
    }
  */

  return 0;
}

