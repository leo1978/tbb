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
#include "vtkStreamTracer.h"
#include "vtkImageGradient.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataWriter.h"

using namespace tbb;
using namespace std;

static const size_t N = 2000;
static const size_t GRAIN = 20;

void CreatePointSource(vtkPoints* in, vtkMultiBlockDataSet* out)
{
  int numBlocks = out->GetNumberOfBlocks();
  int numPoints = in->GetNumberOfPoints();
  int chunkSize = std::max(numPoints/numBlocks,1);

  for(int i=0; i<numBlocks; i++)
    {
    vtkNew<vtkPoints> points_i;
    vtkNew<vtkPolyData> poly_i;
    int chunkStart = i*chunkSize;
    int chunkEnd = min((i+1)*chunkSize,numPoints);
    for(int j=chunkStart; j<chunkEnd; j++)
      {
      points_i->InsertNextPoint(in->GetPoint(j));
      }
    poly_i->SetPoints(points_i.GetPointer());
    out->SetBlock(i, poly_i.GetPointer());
    }
}

vtkSmartPointer<vtkPolyData> MergeStreamLines(vtkMultiBlockDataSet* in)
{
  vtkSmartPointer<vtkAppendPolyData> appender = vtkSmartPointer<vtkAppendPolyData>::New();

  int numBlocks = in->GetNumberOfBlocks();
  for (int idx = 0; idx < numBlocks; ++ idx )
    {
    vtkDataObject* block = in->GetBlock(idx);
    if (block)
      {
      appender->AddInputDataObject(block);
      }
    }
  appender->Update();
  vtkSmartPointer<vtkPolyData> out = vtkPolyData::SafeDownCast(appender->GetOutputDataObject(0));
  return out;
}

int main(int argc, char* argv[])
{
  vtkNew<vtkCompositeDataPipeline> cexec;

  int numRuns = 1;
  int numThreads=0;
  bool noTbb = false;
  int extent = 10;
  int numTraces = 100;
  int numChunks = 0;
  bool dumpOutput = false;
  string filename;

  for(int argi=1; argi<argc; argi++)
    {
    if(string(argv[argi])=="--input")
      {
      filename = string(argv[++argi]);
      }
    else if(string(argv[argi])=="--numThreads")
      {
      numThreads=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--numTraces")
      {
      numTraces=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--noTbb")
      {
      noTbb = true;
      }
    else if(string(argv[argi])=="--dumpOutput")
      {
      dumpOutput = true;
      }
    else if(string(argv[argi])=="--extent")
      {
      extent=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--numRuns")
      {
      numRuns=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--numChunks")
      {
      numChunks=atoi(argv[++argi]);
      }
    else
      {
      cout<<"Say shat?"<<endl;
      return 1;
      }
    }

  bool useThreadedComposite = numThreads>0;

  cout<<"numThreads: "<<numThreads<<endl;
  cout<<numTraces<<" traces"<<endl;
  cout<<"Use "<<(noTbb? "vtk":"tbb")<<" for multithreading"<<endl;
  cout<<"extent: "<<"["<<-extent<<","<<extent<<"]^3"<<endl;
  cout<<(useThreadedComposite? "Use threaded composite executive." : "Use the default executive.")<<endl;
  if(filename!="")
    {
    cout<<filename<<endl;
    }

  if(useThreadedComposite)
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

  if(numThreads>0)
    {
    vtkMultiThreader::SetGlobalDefaultNumberOfThreads(numThreads);
    }

  vtkNew<vtkRTAnalyticSource> source;
  source->SetWholeExtent(-extent,extent,-extent,extent,-extent,extent);

  vtkNew<vtkImageGradient> gradient;
  gradient->SetDimensionality(3);
  gradient->SetInputConnection(source->GetOutputPort());
  gradient->Update();

  vtkNew<vtkPolyData> seeds;
  vtkNew<vtkPoints> seedPoints;

  double d = extent/sqrt((float)numTraces);
  for(double x = -extent/2.0; x< extent/2.0; x+=d)
    {
    for(double y = -extent/2.0; y< extent/2.0; y+=d)
      {
      seedPoints->InsertNextPoint(x,y, 0.0);
      }
    }
  seeds->SetPoints(seedPoints.GetPointer());
  cout<<seeds->GetNumberOfPoints()<<" seeds"<<endl;
  if(numChunks==0)
    {
    numChunks = seeds->GetNumberOfPoints();
    }
  cout<<numChunks<<" chunks"<<endl;
  vtkNew<vtkTimerLog> timer;
  vtkSmartPointer<vtkDataObject> output;
  timer->StartTimer();

  vtkNew<vtkStreamTracer> tracer;
  tracer->SetInputArrayToProcess(0, 0, 0, 0, "RTDataGradient");

  if(!useThreadedComposite)
    {
    tracer->SetSourceData(seeds.GetPointer());
    }
  else
    {
    vtkNew<vtkMultiBlockDataSet> mbSeeds;
    mbSeeds->SetNumberOfBlocks(numChunks);
    CreatePointSource(seedPoints.GetPointer(),mbSeeds.GetPointer());
    tracer->SetInputData(1,mbSeeds.GetPointer());
    }
  tracer->SetInputConnection(gradient->GetOutputPort());
  tracer->SetMaximumPropagation(2000.0);
  tracer->SetMaximumNumberOfSteps(4000);
  tracer->SetIntegrationDirectionToBoth();

  for(int i=1; i<=numRuns; i++)
    {
    if(i==2)
      {
      tracer->Modified();
      }
    tracer->Update();
    if(useThreadedComposite)
      {
      output = MergeStreamLines(vtkMultiBlockDataSet::SafeDownCast(tracer->GetOutputDataObject(0)));
      }
    else
      {
      output = tracer->GetOutputDataObject(0);
      }
    }
  timer->StopTimer();
  std::cout << "time: " << timer->GetElapsedTime() <<" seconds"<< std::endl;
  int numPoints(0);
  int numCells(0);

  if(!vtkCompositeDataSet::SafeDownCast(output))
    {
    cout<<"Report non-composite"<<endl;
    vtkPolyData* out = vtkPolyData::SafeDownCast(output);
    numPoints = out->GetNumberOfPoints();
    numCells = out->GetNumberOfCells();
    if(dumpOutput)
      {
      vtkNew<vtkPolyDataWriter> writer;
      writer->SetInputData(out);
      writer->SetFileName("./stream-traces.vtk");
      writer->Update();
      }
    }
  else
    {
    cout<<"Report composite"<<endl;
    vtkCompositeDataSet* compositeOutput = vtkCompositeDataSet::SafeDownCast(output);
    vtkSmartPointer<vtkCompositeDataIterator> iter;
    iter.TakeReference(compositeOutput->NewIterator());
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
      {
      vtkPolyData* out = vtkPolyData::SafeDownCast(iter->GetCurrentDataObject());
      numPoints+= out->GetNumberOfPoints();
      numCells+= out->GetNumberOfCells();
      }
    }

  cout<<numPoints<<" points, "<<numCells<<" cells "<<endl;

  vtkAlgorithm::SetDefaultExecutivePrototype(NULL);

  return 0;
}

