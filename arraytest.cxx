#include "vtkCellArray.h"
#include <vtkTimerLog.h>
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkCompositeDataIterator.h"
#include "vtkPointData.h"
#include "vtkPoints.h"

#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkBoundingBox.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>

using namespace std;

struct vec3
{
  float x[3];
  vec3(float a, float b, float c)
  {
    x[0] = a;
    x[1] = b;
    x[2] = c;
  }
};

int main(int argc, char* argv[])
{
  int n(1000000);
  int times(10);
  bool useStl=false;
  bool useHelper = false;
  bool useC=false;
  bool check = false;
  bool useLazyCopy =false;
  string op = "insert";
  for(int argi=1; argi<argc; argi++)
    {
    if(string(argv[argi])=="--op")
      {
      op = argv[++argi];
      }
    if(string(argv[argi])=="--times")
      {
      times=atoi(argv[++argi]);
      }
    if(string(argv[argi])=="--stl")
      {
      useStl=true;
      }
    if(string(argv[argi])=="--helper")
      {
      useHelper=true;
      }
    if(string(argv[argi])=="--lazy")
      {
      useLazyCopy=true;
      }
    if(string(argv[argi])=="--check")
      {
      check=true;
      }
    if(string(argv[argi])=="--c")
      {
      useC=true;
      }
    if(string(argv[argi])=="--n")
      {
      n=atoi(argv[++argi]);
      }
    }
  vtkPointData* pd = vtkPointData::New();
  std::vector<std::vector<float> > pd1(3);
  if(op==string("copy"))
    {
      {
      vtkFloatArray* A = vtkFloatArray::New();
      vtkFloatArray* B = vtkFloatArray::New();
      vtkFloatArray* C = vtkFloatArray::New();
      pd->AddArray(A); A->Delete();
      pd->AddArray(B); B->Delete();
      pd->AddArray(C); C->Delete();
      for(int i=0; i<n; i++)
        {
        A->InsertValue(i, (float)i);
        B->InsertValue(i, (float)i);
        C->InsertValue(i, (float)i);
        }
      }
      {
      std::vector<float>& A(pd1[0]);
      std::vector<float>& B(pd1[1]);
      std::vector<float>& C(pd1[2]);
      for(int i=0; i<n; i++)
        {
        A.push_back((float)i);
        B.push_back((float)i);
        C.push_back((float)i);
        }
      }
    }

  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  for(int round=0; round<times; round++)
    {
    if(op=="insert")
      {
      if(useStl)
        {
        std::vector<vec3> points;
        for(int i=0; i<n; i++)
          {
          vec3 a(0,1,2);
          points.push_back(a);
          }
        }
      else if (useC)
        {
        float* points = new float[3*n];
        for(int i=0; i<n; i++)
          {
          float x[3] = {0,1,2};
          points[3*i] = x[0];
          points[3*i+1] = x[1];
          points[3*i+2] = x[2];
          }
        delete points;
        }
      else
        {
        vtkPoints* points = vtkPoints::New();
        for(int i=0; i<n; i++)
          {
          float x[3]={0,1,2};
          points->InsertNextPoint(x);
          }
        assert(points->GetNumberOfPoints()==n);
        points->Delete();
        }
      }
    else if(op=="copy")
      {
      vtkPointData* outPd = vtkPointData::New();
      if(useStl)
        {
        std::vector<std::vector<float> > outPd(3);
        for(int i=0; i<3; i++)
          {
          std::vector<float>& src = pd1[i];
          std::vector<float>& dst = outPd[i];
          for(int j=0; j<n; j++)
            {
            dst.push_back(src[j]);
            }
          }
        }
      else if (useC)
        {
        float* outPD[3];
        for(int i=0; i<3; i++)
          {
          outPD[i] = new float[n];
          memcpy(outPD[i], &pd1[i][0], sizeof(float)*n);
          }
        for(int i=0; i<3; i++)
          {
          delete[] outPD[i];
          }
        }
      else if (useHelper)
        {
        outPd->CopyAllocate(pd); //MUST DO THIS
        CopyHelper copy(pd,outPd);
        for(int i=0; i<n; i++)
          {
          copy.Copy(i,i);
          }
        }
      else if (useLazyCopy)
        {
        outPd->CopyAllocate(pd); //MUST DO THIS
        LazyCopy copy(pd,outPd);
        for(int i=0; i<n; i++)
          {
          copy.Copy(i,i);
          }
        copy.Flush();
        }
      else
        {
        outPd->CopyAllocate(pd);
        for(int i=0; i<n; i++)
          {
          outPd->CopyData(pd, i, i);
          }
        }
      if(check)
        {
        for(int i=0; i<n; i++)
          {
          for(int j=0; j<3; j++)
            {
            float value = vtkFloatArray::SafeDownCast(outPd->GetArray(j))->GetValue(i);
            if(!(vtkFloatArray::SafeDownCast(outPd->GetArray(j))->GetValue(i)== (float)i))
              {
              cout<<"Results are wrong."<<endl;
              return 1;
              }
            }
          }
        }
      outPd->Delete();
      }
    }

  timer->StopTimer();
  cout<<"Time: "<<timer->GetElapsedTime()<<" seconds"<<endl;

  pd->Delete();
  return 0;
}

