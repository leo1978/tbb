#include <iostream>
#include <string>
#include <algorithm>
#include <tbb/tbb.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "vtkCellDataToPointData.h"
#include "vtkCellArray.h"

#include <vtkSynchronizedTemplates3D.h>
#include <vtkSynchronizedTemplatesCutter3D.h>

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
#include "vtkCompositeDataIterator.h"
#include "vtkPointData.h"
#include "vtkExtentTranslator.h"
#include "vtkPolyDataWriter.h"
#include "vtkAppendPolyData.h"
#include "vtkCleanPolyData.h"
#include "vtkMergePoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGeometryFilter.h"
#include "vtkPointDataToCellData.h"

#include "vtkContourFilter.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkContourGrid.h"
#include "vtkThreshold.h"
#include "vtkPointSetNormalEstimation.h"
#include "vtkCurvatures.h"
#include "vtkMarchingCubes.h"
#include "vtkFloatArray.h"

#include <algorithm>

using namespace tbb;
using namespace std;

static const size_t N = 2000;
static const size_t GRAIN = 20;


const float supernova_isovalue=0.07;

void ReadSuperNova(vtkImageData* image)
{
  const int GRID_SIZE = 432;
  try
    {
    //read the super nova data from file -> buffer
    float* buffer = new float[GRID_SIZE*GRID_SIZE*GRID_SIZE];
    FILE *fd = fopen("/home/leo/MarchingCubesComparison/SuperNova/normal_1349.dat", "rb");
    assert(fd != NULL);
    size_t numPixels = fread(buffer, sizeof(float), GRID_SIZE*GRID_SIZE*GRID_SIZE, fd);
    assert(ferror(fd) == 0);
    assert(numPixels == GRID_SIZE*GRID_SIZE*GRID_SIZE);
    fclose(fd);
    cout<<"range: "<<*(std::min_element(buffer,buffer+numPixels))<<" "
        <<*(std::max_element(buffer,buffer+numPixels))<<endl;

    //buffer -> image
    image->SetOrigin(0.0, 0.0, 0.0);
    image->SetSpacing(1.0, 1.0, 1.0);
    image->SetExtent(0, GRID_SIZE-1,0, GRID_SIZE-1,0, GRID_SIZE-1);

    vtkNew<vtkFloatArray> vtkElevationPoints;
    vtkElevationPoints->SetName("Elevation");
    vtkElevationPoints->SetVoidArray(&buffer[0],numPixels,1);

    image->GetPointData()->AddArray(vtkElevationPoints.GetPointer());
    image->GetPointData()->SetActiveScalars("Elevation");
    }
  catch (std::exception &error)
    {
    std::cout << "Caught standard exception: " << error.what() << std::endl;
    }

}

void CreatePointSource(vtkPoints* in, vtkMultiBlockDataSet* out)
{
  int numPoints = in->GetNumberOfPoints();
  out->SetNumberOfBlocks( numPoints );

  for(int i=0; i<numPoints; i++)
    {
    vtkNew<vtkPoints> points_i;
    vtkNew<vtkPolyData> poly_i;
    points_i->InsertNextPoint(in->GetPoint(i));

    poly_i->SetPoints(points_i.GetPointer());
    out->SetBlock(i, poly_i.GetPointer());
    }
}


//----------------------------------------------------------------------------
// returns the next pointer in dest
vtkIdType *AppendCells(vtkIdType *pDest, vtkCellArray *src, vtkIdType* vertexMap)
{
  vtkIdType *pSrc, *end, *pNum;

  if (src == NULL)
    {
    return pDest;
    }

  pSrc = src->GetPointer();
  end = pSrc + src->GetNumberOfConnectivityEntries();
  pNum = pSrc;

  while (pSrc < end)
    {
    if (pSrc == pNum)
      {
      // move cell pointer to next cell
      pNum += 1+*pSrc;
      // copy the number of cells
      *pDest++ = *pSrc++;
      }
    else
      {
      // offset the point index
      *pDest++ = vertexMap[*pSrc++];
      }
    }

  return pDest;
}

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkBoundingBox.h"
inline void Copy(double in[3], double out[3])
{
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
}

void DuplicateDataObject(vtkInformation* info, int n, vtkInformationVector* infos)
{
  vtkDataObject* data = info->Get(vtkDataObject::DATA_OBJECT());

  infos->SetNumberOfInformationObjects(n);
  for(int i=0; i<n; i++)
   {
   //Allocate a new information object
   vtkNew<vtkInformation> info_i;
   info_i->Copy(info,1);
   infos->SetInformationObject(i, info_i.GetPointer());

   //Fill it out
   vtkDataObject* data_i = data->NewInstance();
   data_i->ShallowCopy(data);
   info_i->Set(vtkDataObject::DATA_OBJECT(), data_i);
   data_i->FastDelete();
   }
}


struct TaskChunks
{
  int Port;
  int Conn;
  vtkSmartPointer<vtkInformationVector> Input;
  vtkSmartPointer<vtkInformationVector> Output;
  TaskChunks():Port(-1),Conn(-1)
  {
    this->Input = vtkSmartPointer<vtkInformationVector>::New();
    this->Output = vtkSmartPointer<vtkInformationVector>::New();
  }
  int Size()
  {
    assert(this->Input->GetNumberOfInformationObjects()==this->Output->GetNumberOfInformationObjects());
    return this->Input->GetNumberOfInformationObjects();
  }
};

namespace
{
  static vtkInformationVector** Clone(vtkInformationVector** src, int n)
  {
    vtkInformationVector** dst = new vtkInformationVector*[n];
    for(int i=0; i<n; ++i)
      {
      dst[i] = vtkInformationVector::New();
      dst[i]->Copy(src[i],1);
      }
    return dst;
  }
  static void DeleteAll(vtkInformationVector** dst, int n)
  {
    for(int i=0; i<n; ++i)
      {
      dst[i]->Delete();
      }
    delete []dst;
  }
};

class FilterTask
{
public:
  FilterTask(vtkAlgorithm* filter, TaskChunks& chunks):Filter(filter),Chunks(chunks)
  {
    this->Input = this->Filter->GetExecutive()->GetInputInformation();
    this->Output = this->Filter->GetExecutive()->GetOutputInformation();
  }
  void operator() ( const blocked_range<int>& r ) const
  {
    vtkSmartPointer<vtkInformation> requestData = vtkSmartPointer<vtkInformation>::New();
    requestData->Set(vtkDemandDrivenPipeline::REQUEST_DATA());

    int numPorts = this->Filter->GetNumberOfInputPorts();
    vtkInformationVector** input = Clone(this->Input, numPorts );
    vtkInformationVector* output = vtkInformationVector::New();
    output->Copy(this->Output,1);

    for ( int i = r.begin(); i != r.end(); ++i )
      {
      vtkInformation* pInput = this->Chunks.Input->GetInformationObject(i);
      vtkInformation* pOutput = this->Chunks.Output->GetInformationObject(i);
      input[this->Chunks.Port]->SetInformationObject(this->Chunks.Conn, pInput);
      output->SetInformationObject(0, pOutput);
      this->Filter->ProcessRequest(requestData, input, output);
      }

    DeleteAll(input,numPorts);
    output->Delete();
   }

private:
  vtkAlgorithm* Filter;
  TaskChunks& Chunks;
  vtkInformationVector** Input;
  vtkInformationVector* Output;
};

vtkSmartPointer<vtkAlgorithm> CreateTestFilter(const string& filterName)
{
  if(filterName ==string("sync"))
    {
    vtkNew<vtkSynchronizedTemplates3D> contour;
    double contourValue(157.0);
    contour->SetValue(0,contourValue);
    return contour.GetPointer();
    }
  else if(filterName ==string("contour"))
    {
    vtkNew<vtkContourFilter> contour;
    double contourValue(157.0);
    contour->SetValue(0,contourValue);
    return contour.GetPointer();
    }
  else if(filterName ==string("mc"))
    {
    vtkNew<vtkMarchingCubes> contour;
    double contourValue(157.0);
    contour->SetValue(0,contourValue);
    return contour.GetPointer();
    }
  else if(filterName ==string("cgrid"))
    {
    vtkNew<vtkContourGrid> contour;
    return contour.GetPointer();
    }
  else if(filterName ==string("thresh"))
    {
    vtkNew<vtkThreshold> thresh;
    return thresh.GetPointer();
    }

  else if(filterName ==string("cutter"))
    {
    vtkNew<vtkSynchronizedTemplatesCutter3D> cutter;
    vtkNew<vtkPlane> plane;
    plane->SetOrigin(0,0,0);
    plane->SetNormal(1,1,1);
    cutter->SetCutFunction(plane.GetPointer());
    return cutter.GetPointer();
    }
  else if(filterName ==string("p2c"))
    {
    vtkNew<vtkPointDataToCellData> filter;
    return filter.GetPointer();
    }
  else if(filterName ==string("normal"))
    {
    vtkNew<vtkPointSetNormalEstimation> filter;
    return filter.GetPointer();
    }
  else if(filterName ==string("surf"))
    {
    vtkNew<vtkGeometryFilter> filter;
    return filter.GetPointer();
    }
  else
    {
    cout<<"Does not understand what a "<<filterName<<" is"<<endl;
    return NULL;
    }
}

int main(int argc, char* argv[])
{
  int numThreads=0;
  int extent = 10;
  bool dumpOutput = false;
  int numChunks(64);
  bool noTbb(false);
  bool superNova(false);
  bool noMerge(false);
  string filterName = "sync";

  for(int argi=1; argi<argc; argi++)
    {
    if(string(argv[argi])=="--numThreads")
      {
      numThreads=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--dumpOutput")
      {
      dumpOutput = true;
      }
    else if(string(argv[argi])=="--extent")
      {
      extent=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--numChunks")
      {
      numChunks=atoi(argv[++argi]);
      }
    else if(string(argv[argi])=="--noTbb")
      {
      noTbb = true;
      }
    else if(string(argv[argi])=="--noMerge")
      {
      noMerge = true;
      }
    else if(string(argv[argi])=="--filter")
      {
      filterName = string(argv[++argi]);
      }
    else if(string(argv[argi])=="--superNova")
      {
      superNova = true;
      }
    else
      {
      cout<<"Say shat?"<<endl;
      return 1;
      }
    }
  cout<<numThreads<<" threads"<<endl;
  cout<<numChunks<<" pieces"<<endl;
  vtkThreadedCompositeDataPipeline::NumChunks = numChunks;

  //chosen image source -> image
  vtkSmartPointer<vtkImageData> image;
  if(superNova)
    {
    image = vtkSmartPointer<vtkImageData>::New();
    ReadSuperNova(image);
    }
  else
    {
    vtkNew<vtkRTAnalyticSource> source;
    source->SetWholeExtent(-extent,extent,-extent,extent,-extent,extent);
    source->Update();
    image = vtkImageData::SafeDownCast(source->GetOutputDataObject(0));
    assert(image->GetPointData()->GetAbstractArray("RTData"));
    }

  vtkNew<vtkFloatArray> arr;
  arr->SetNumberOfTuples(image->GetNumberOfCells());
  arr->SetNumberOfComponents(1);
  arr->FillComponent(0,0);
  image->GetCellData()->AddArray(arr.GetPointer());

  vtkSmartPointer<vtkAlgorithm> filter = CreateTestFilter(filterName);

  filter->SetInputDataObject(0,image);

  if(filterName=="sync")
    {
    vtkSynchronizedTemplates3D* contour = vtkSynchronizedTemplates3D::SafeDownCast(filter);
    assert(contour);
    contour->SetNumberOfContours(1);
    contour->ComputeScalarsOn();
    double contourValue = superNova? supernova_isovalue : 157.0;
    cout<<"Set contour value: "<<contourValue<<endl;
    contour->SetValue(0,contourValue);
    }

  if(superNova)
    {
    filter->SetInputArrayToProcess(0, 0, 0, 0, "Elevation");
    }
  else
    {
    filter->SetInputArrayToProcess(0, 0, 0, 0, "RTData");
    }


  int numPoints(0);
  int numCells(0);
  int numCellData(0);
  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();

  vtkSmartPointer<vtkPolyData> output;
  if(numThreads==0)
    {
    filter->Update();
    vtkDataObject* out = filter->GetOutputDataObject(0);
    output = vtkPolyData::SafeDownCast(out);
    }
  else
    {

    //new
    filter->UpdateInformation();
    filter->PropagateUpdateExtent();

    TaskChunks chunks;
    chunks.Port = filter->GetParallelInputPort();
    chunks.Conn = 0;
    vtkSmartPointer<vtkInformation> requestChunks = vtkSmartPointer<vtkInformation>::New();
    requestChunks->Set(vtkThreadedCompositeDataPipeline::REQUEST_DIVIDE());
    filter->ProcessRequest(requestChunks, filter->GetExecutive()->GetInputInformation(),chunks.Input);
    int numChunks = chunks.Input->GetNumberOfInformationObjects();
    DuplicateDataObject(filter->GetOutputInformation(0), numChunks,chunks.Output);
    assert(chunks.Size()==numChunks);

    tbb::blocked_range<int> range(0, chunks.Size());
    FilterTask task(filter.GetPointer(),chunks);

    if(noTbb)
      {
      vtkMultiThreader::SetGlobalDefaultNumberOfThreads(numThreads);
      vtkParallelFor(range,task);
      }
    else
      {
      tbb::task_scheduler_init init(numThreads);
      parallel_for(range,task);
      }

    vtkNew<vtkTimerLog> timer1;
    timer1->StartTimer();
    if(!noMerge)
      {
      vtkSmartPointer<vtkInformation> request = vtkSmartPointer<vtkInformation>::New();
      request->Set(vtkThreadedCompositeDataPipeline::REQUEST_MERGE());

      vtkInformationVector* infoVec[1] = {chunks.Output};
      filter->ProcessRequest(request, &infoVec[0], filter->GetExecutive()->GetOutputInformation());
      output = vtkPolyData::GetData(filter->GetOutputInformation(0));

      timer1->StopTimer();
      }
    else
      {
      // vtkSmartPointer<vtkAppendPolyData> appender = vtkSmartPointer<vtkAppendPolyData>::New();
      // for (int idx = 0; idx < numChunks; ++ idx )
      //   {
      //   vtkInformation* outInfo = chunks.Output->GetInformationObject(idx);
      //   vtkPolyData* outi = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
      //   appender->AddInputDataObject(outi);
      //   }
      // appender->Update();
      // output = appender->GetOutput();
      }
    cout<<"Merge took "<<timer1->GetElapsedTime()<<" seconds"<<endl;
    }

  timer->StopTimer();

  cout<<"Report output"<<endl;
  if(vtkUnstructuredGrid::SafeDownCast(filter->GetOutputDataObject(0)))
    {
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(filter->GetOutputDataObject(0));
    numPoints = grid->GetNumberOfPoints();
    numCells = grid->GetNumberOfCells();
    numCellData = grid->GetCellData()->GetNumberOfTuples();
    }
  else if(vtkImageData::SafeDownCast(filter->GetOutputDataObject(0)))
    {
    vtkImageData* image = vtkImageData::SafeDownCast(filter->GetOutputDataObject(0));
    numPoints = image->GetNumberOfPoints();
    numCells = image->GetNumberOfCells();
    numCellData = image->GetCellData()->GetNumberOfTuples();
    }
  else
    {
    numPoints = output? output->GetNumberOfPoints() : 0;
    numCells = output? output->GetNumberOfCells() : 0;
    numCellData = output? output->GetCellData()->GetNumberOfTuples():0;
    }

  std::cout << "time: " << timer->GetElapsedTime() <<" seconds"<< std::endl;

  cout<<numPoints<<" points, "<<numCells<<" cells "<<numCellData<<" cell data values"<<endl;

  if(dumpOutput)
    {
    char fname[] = "./contour.vtk";
    cout<<"Write to :"<<fname<<endl;
    vtkNew<vtkPolyDataWriter> writer;
    writer->SetInputData(output);
    writer->SetFileName(fname);
    writer->Update();
    }

  return 0;
}

