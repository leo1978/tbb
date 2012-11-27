#include <iostream>
#include <string>
#include <algorithm>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "vtkCellDataToPointData.h"

#include <vtkSynchronizedTemplates3D.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkRTAnalyticSource.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include "vtkCompositeDataPipeline.h"
#include "vtkThreadedCompositeDataPipeline.h"
#include <vtkTimerLog.h>
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkParallelFor.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMutexLock.h"
#include "vtkAMREnzoReader.h"
#include "vtkAMRFlashReader.h"
#include "vtkCompositeDataIterator.h"
using namespace tbb;
using namespace std;

static const size_t N = 2000;
static const size_t GRAIN = 20;

int main(int argc, char* argv[])
{
  vtkNew<vtkCompositeDataPipeline> cexec;

  bool bigData = false;
  bool threadedCompositePipeline = true;
  int numThreads=0;
  bool noTbb = false;
  string filename;

  for(int argi=1; argi<argc; argi++)
    {
    if(string(argv[argi])=="--bigdata")
      {
      bigData = true;
      }
    else if(string(argv[argi])=="--input")
      {
      filename = string(argv[++argi]);
      }
    else if(string(argv[argi])=="--compositePipeline")
      {
      threadedCompositePipeline=false;
      }
    else if(string(argv[argi])=="--numThreads")
      {
      numThreads=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--noTbb")
      {
      noTbb = true;
      }
    else
      {
      cout<<"Say shat?"<<endl;
      return 1;
      }
    }

  //not doing multithreading by default
  if(numThreads==0)
    {
    threadedCompositePipeline = false;
    }

  string fname;
  string scalar_name;
  vtkSmartPointer<vtkAMRBaseReader> reader;
  double contourValue(0);
  if(bigData)
    {
    fname = "/home/leo/pvdata/smooth.flash";
    scalar_name = "dens";
    reader = vtkSmartPointer<vtkAMRFlashReader>::New();
    contourValue = 3000.0;
    }
  else
    {
    scalar_name = "Density";
    fname = "/home/leo/VTKData/Data/AMR/Enzo/DD0010/moving7_0010.hierarchy";
    reader = vtkSmartPointer<vtkAMREnzoReader>::New();
    vtkAMREnzoReader::SafeDownCast(reader)->ConvertToCGSOff();
    contourValue = 30;
    }


  cout<<"data set: "<<fname<<", "<<scalar_name<<endl;
  cout<<"threadedCompositePipeline "<<threadedCompositePipeline<<endl;
  cout<<"numThreads: "<<numThreads<<endl;
  cout<<"noTbb: "<<noTbb<<endl;
  if(filename!="")
    {
    cout<<filename<<endl;
    }

  if(numThreads>0)
    {
    vtkMultiThreader::SetGlobalDefaultNumberOfThreads(numThreads);
    }
  reader->SetFileName(fname.c_str());
  reader->SetCellArrayStatus(scalar_name.c_str(),1);
  reader->SetMaxLevel(5);
  reader->Update();
  cout<<"Done reading AMR"<<endl;
  int total_blocks(0);
    {
    vtkCompositeDataSet* output = vtkCompositeDataSet::SafeDownCast(reader->GetOutputDataObject(0));
    vtkSmartPointer<vtkCompositeDataIterator> iter;
    iter.TakeReference(output->NewIterator());
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
      {
      if(iter->GetCurrentDataObject())
        {
        total_blocks++;
        }
      }
    cout<<total_blocks<<" blocks"<<endl;
    }

  if(threadedCompositePipeline)
    {
    vtkNew<vtkThreadedCompositeDataPipeline> exec;
    vtkThreadedCompositeDataPipeline::UseTBB = !noTbb;
    vtkAlgorithm::SetDefaultExecutivePrototype(exec.GetPointer());
    }
  else
    {
    vtkNew<vtkCompositeDataPipeline> exec;
    vtkAlgorithm::SetDefaultExecutivePrototype(exec.GetPointer());
    }

  vtkNew<vtkCellDataToPointData> c2p;
  c2p->SetInputConnection(reader->GetOutputPort());

  vtkNew<vtkSynchronizedTemplates3D> cf;
  cf->SetInputConnection(c2p->GetOutputPort());
  cf->SetInputArrayToProcess(0, 0, 0, 0, scalar_name.c_str());
  cf->SetValue(0,contourValue);

  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  cf->Update();
  timer->StopTimer();
  std::cout << "time: " << timer->GetElapsedTime() <<" seconds"<< std::endl;

  int numPolygons(0);
  vtkCompositeDataSet* output = vtkCompositeDataSet::SafeDownCast(cf->GetOutputDataObject(0));
  vtkSmartPointer<vtkCompositeDataIterator> iter;
  iter.TakeReference(output->NewIterator());
  int numBlocks(0);
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
    vtkPolyData* out = vtkPolyData::SafeDownCast(iter->GetCurrentDataObject());
    numPolygons+= out->GetNumberOfCells();
    numBlocks++;
    }

  cout<<numPolygons<<" polygons in "<<numBlocks<<" blocks"<<endl;

  vtkAlgorithm::SetDefaultExecutivePrototype(NULL);

  return 0;
}
