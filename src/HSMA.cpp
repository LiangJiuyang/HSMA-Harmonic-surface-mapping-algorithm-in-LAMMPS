 /* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "math.h"
#include "fmm.h"
#include<iostream>
#include"mkl.h"
#include<omp.h>
#include<iomanip>

#include "complex.h"

extern "C" {void lfmm3d_t_c_g_(double *eps, int *nsource,double *source, double *charge, int *nt, double *targ, double *pottarg, double *gradtarg, int *ier);}


using namespace LAMMPS_NS; 
using namespace std;

/* ---------------------------------------------------------------------- */
HSMA::HSMA(LAMMPS *lmp) : KSpace(lmp)
{
  maxatom = atom->natoms;
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &RankID);
  Lx = domain->xprd;
  Ly = domain->yprd;
  Lz = domain->zprd * slab_volfactor;
  pi = 3.141592653589793;

  epot = NULL;
  efield = NULL;

  //ewaldflag = 1;
}

/* ---------------------------------------------------------------------- */
void HSMA::settings(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Illegal kspace_style FMM command");
  tolerance = fabs(force->numeric(FLERR,arg[0]));

  Lambda = force->numeric(FLERR, arg[1]);
  p = force->numeric(FLERR, arg[2]);
  Nw = force->numeric(FLERR, arg[3]);
  Fp = force->numeric(FLERR, arg[4]);
  F = force->numeric(FLERR, arg[5]);
  IF_FMM_RightTerm = force->numeric(FLERR, arg[6]);
  IF_FMM_FinalPotential = force->numeric(FLERR, arg[7]);
}

void HSMA::init()
{
  printf("Setting up FMM implemented by Greengard (Release 1.0.0)\n");

  //设置测试点的位置
  //double srun1[Nw][3],srun2[Nw][3];
  //PointSum = srun1;
  //QuizSum = srun2;
  
  PointSum = new double * [Nw];
  QuizSum = new double * [Nw];
  for (int i = 0; i < Nw; i++)
  {
      PointSum[i] = new double[3];
      QuizSum[i] = new double[3];
  }

  R= sqrt((Lx / 2) * (Lx / 2) + (Ly / 2) * (Ly / 2) + (Lz / 2) * (Lz / 2));
  Rs = Lambda* sqrt((Lx / 2) * (Lx / 2) + (Ly / 2) * (Ly / 2) + (Lz / 2) * (Lz / 2));//The radius of the cut-off sphere B
  for (int i = 0; i < Nw; i++)
  {
	  PointSum[i][2] = (2 * (i + 1) - 1) / (Nw + 0.00) - 1;
	  PointSum[i][0] = (sqrt(1 - PointSum[i][2] * PointSum[i][2]) * cos(2 * pi * (i + 1) * 0.618)) * R;
	  PointSum[i][1] = (sqrt(1 - PointSum[i][2] * PointSum[i][2]) * sin(2 * pi * (i + 1) * 0.618)) * R;
	  PointSum[i][2] = PointSum[i][2] * R;


	  if (abs(PointSum[i][0]) >= (Lx / 2))
	  {
		  QuizSum[i][0] = (abs(PointSum[i][0]) - Lx * int((abs(PointSum[i][0]) + Lx / 2) / Lx)) * PointSum[i][0] / abs(PointSum[i][0]);
	  }
	  else
	  {
		  QuizSum[i][0] = PointSum[i][0];
	  }
	  if (abs(PointSum[i][1]) >= (Ly / 2))
	  {
		  QuizSum[i][1] = (abs(PointSum[i][1]) - Ly * int((abs(PointSum[i][1]) + Ly / 2) / Ly)) * PointSum[i][1] / abs(PointSum[i][1]);
	  }
	  else
	  {
		  QuizSum[i][1] = PointSum[i][1];
	  }
	  if (abs(PointSum[i][2]) >= (Lz / 2))
	  {
		  QuizSum[i][2] = (abs(PointSum[i][2]) - Lz * int((abs(PointSum[i][2]) + Lz / 2) / Lz)) * PointSum[i][2] / abs(PointSum[i][2]);
	  }
	  else
	  {
		  QuizSum[i][2] = PointSum[i][2];
	  }
  }

  //计算测试点的局部展开基
  A = new double* [Nw];
  PointSumMultipleExpansionMatrix = new double[p * p];
  QuizSumMultipleExpansionMatrix = new double[p * p];
  for (int i = 0; i < Nw; i++)
  {
	  A[i] = new double[p*p-1];
  }
  for (int i = 0; i < Nw; i++)
  {
	  CalculateMultipoleExpansion(PointSumMultipleExpansionMatrix, p, PointSum[i][0], PointSum[i][1], PointSum[i][2]);
	  CalculateMultipoleExpansion(QuizSumMultipleExpansionMatrix, p, QuizSum[i][0], QuizSum[i][1], QuizSum[i][2]);
	  for (int j = 0; j < p * p - 1; j++)
	  {
		  A[i][j] = PointSumMultipleExpansionMatrix[j + 1] - QuizSumMultipleExpansionMatrix[j + 1];
	  }
  }

  //Construct Fibonacci Points And There Local/Multipole Expansion.  Fibonacci[Np][4] storage the positions and the weights.
  Np = int(2 * F + 2.0);
  Fibonacci = new double* [Np];
  for (int i = 0; i < Np; i++)
  {
	  Fibonacci[i] = new double[4];
  }
  double Fibonacci_New[Np][4];
  SetFibonacci(Fibonacci_New, F, Fp, Np, Rs, pi);
  
  for (int i = 0; i < Np; i++)
	  for (int j = 0; j < 4; j++)
		  Fibonacci[i][j] = Fibonacci_New[i][j];

  QRD = new double* [Np]; QLocalRD= new double* [Np];
  for (int i = 0; i < Np; i++)
  {
	  QRD[i] = new double[p*p];
	  QLocalRD[i]= new double[p * p];
  }
  #pragma omp parallel
  {
	  double MulQ[p * p], MulLocalQ[p * p];
      #pragma omp for schedule(static) private(MulQ,MulLocalQ)
	  for (int i = 0; i < Np; i++)
	  {
		  CalculateRDMultipoleExpansion(MulQ, p, Fibonacci[i][0], Fibonacci[i][1], Fibonacci[i][2]);
		  CalculateLocalRDMultipoleExpansion(MulLocalQ, p, Fibonacci[i][0], Fibonacci[i][1], Fibonacci[i][2], Rs);
		  for (int j = 0; j < p * p; j++)
		  {
			  QRD[i][j] = MulQ[j];
			  QLocalRD[i][j] = MulLocalQ[j];
		  }
	  }
  }
  cout << Lx << "=Lx   " << Ly << "=Ly    " << Lz << "=Lz   " << Lambda << "=Lambda    " << p << "=p   " << Nw << "=Nw    "  << Fp << "=Fp   " << F << "=F   " << IF_FMM_RightTerm << "=IF_FMM_RightTerm   " << IF_FMM_FinalPotential << "=IF_FMM_FinalPotential   " << endl;
}


void HSMA::compute(int eflag, int vflag)
{
	// set energy/virial flags
	ev_init(eflag, vflag);
	// if atom count has changed, update qsum and qsqsum
	if (atom->natoms != natoms_original) {
		qsum_qsq();
		natoms_original = atom->natoms;
	}
	// return if there are no charges
	if (qsqsum == 0.0) return;

  //设置粒子数 坐标 带电量 力的接口	
  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  qqrd2e = force->qqrd2e;
  double **f = atom->f;
  double boxlo[3] = { domain->boxlo[0] ,domain->boxlo[1],domain->boxlo[2]};
  double boxhi[3] = { domain->boxhi[0] ,domain->boxhi[1],domain->boxhi[2]};



  if (me == 0) {
	  //cout << "The messages of box are " << domain->boxlo[0] << "   " << domain->boxhi[0] << endl;
	  //cout << "The messages of box are " << domain->boxlo[1] << "   " << domain->boxhi[1] << endl;
	  //cout << "The messages of box are " << domain->boxlo[2] << "   " << domain->boxhi[2] << endl;
	  //cout << C[0] << "   " << C[1] << "    " << C[2] << endl;
	  //cout << PointSum[10][0] << "   " << PointSum[10][1] << "    " << PointSum[10][2] << endl;
	  //cout << QuizSum[10][0] << "   " << QuizSum[10][1] << "    " << QuizSum[10][2] << endl;
	  //cout << Near[0] << "   " << Near[1] << "    " << Near[2] << endl;
	  //cout << Near_All[0] << "   " << Near_All[1] << "    " << Near_All[2] << endl;
	  //cout << x[0][0] << "   " << x[0][1] << "   " << x[0][2] << endl;
	  //cout << f[0][0] << "   " << f[0][1] << "   " << f[0][2] << endl;
  }

  if (RankID == 1)
  {
	  double time;
	  time = MPI_Wtime();

	  double X[nlocal][3], Q[nlocal], Force[nlocal][3], Pot[nlocal];
	  for (int i = 0; i < nlocal; i++)
	  {
		  X[i][0] = x[i][0] - (boxhi[0] - Lx / 2);
		  X[i][1] = x[i][1] - (boxhi[1] - Ly / 2);
		  X[i][2] = x[i][2] - (boxhi[2] - Lz / 2);
		  Q[i] = q[i];
		  Force[i][0] = 0.00;
		  Force[i][1] = 0.00;
		  Force[i][2] = 0.00;
		  Pot[i] = 0.00;
	  }
	  AdjustParticle_Double(X, nlocal, Lx, Ly, Lz);

	  //time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for Final is " << time << endl;
	  //time = MPI_Wtime();

	  //Find the image charge(这里每个核只存储局部的镜像)
	  int lx = ceil((Rs - Lx / 2) / Lx), ly = ceil((Rs - Ly / 2) / Ly), lz = ceil((Rs - Lz / 2) / Lz);
	  int TotalNumber = ceil(nlocal * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);//减少存储量消耗
	  double ImageCharge[TotalNumber][4];
	  int ImageNumber;
	  SetImageCharge(ImageCharge, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), X, Q, nlocal, Rs, Lx, Ly, Lz, lx, ly, lz);

	  //time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for Final is " << time << endl;
	  //time = MPI_Wtime();

	  //Calculate near field potential which contains all the contributions in the cut-off sphere B at the monitoring points
	  double Near[Nw];
	  CalculateNearFieldAndZD_Single(Near, ImageCharge, ImageNumber, Nw, IF_FMM_RightTerm, Lx, Ly, Lz, PointSum, QuizSum, X, Force, Pot, nlocal, Q,tolerance);

	  //time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for Final is " << time << endl;
	  //time = MPI_Wtime();

	  //Solve The Least Square Problem
	  double C[p * p];
	  SolveLeastSquareProblem(C, A, Near, p, Nw);

	  //time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for Final is " << time << endl;
	  //time = MPI_Wtime();

	  //Compute final force, potential and energy
	  double Energy_HSMA;
	  Energy_HSMA = FinalCalculateEnergyAndForce_Single(Force, Pot, X, Q, nlocal, ImageCharge, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, C, p, Fp, F, Rs, pi, IF_FMM_FinalPotential,tolerance);

	  time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for Final is " << time << endl;

	  scale = 1.0;
	  const double qscale = qqrd2e * scale;

	  //Assemble results
	  for (int i = 0; i < nlocal; i++)
	  {
		  f[i][0] += Force[i][0] * qscale;
		  f[i][1] += Force[i][1] * qscale;
		  f[i][2] += Force[i][2] * qscale;
	  }

	  if (eflag_global) {
		  energy += Energy_HSMA * qscale / 2;
		  //energy += Energy_HSMA * qscale; 
	  }

	  if (evflag_atom) {
		  if (eflag_atom) {
			  for (int i = 0; i < nlocal; i++) {
				  eatom[i] += Pot[i];
				  eatom[i] *= qscale;
			  }
		  }

		  if (vflag_atom)
			  for (int i = 0; i < nlocal; i++)
				  for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
	  }

	  if (me == 0) {
		  //cout << "The forces of particle 0 are" << f[0][0] << "   " << f[0][1] << "   " << f[0][2] << endl;
		  //cout << "The Pot of particle 0 are" << Pot[0] << endl;
		  //cout << "The Energy is " << energy << "      " << Energy_HSMA / maxatom << endl;
	  }

  }
  else if (RankID > 1) {

	  double AllSource[maxatom][3], AllQ[maxatom];

	  //收集全部粒子的信息
	  int nlocal_All[RankID], nlocal_All_Q[RankID];
	  MPI_Allgather(&nlocal, 1, MPI_INT, nlocal_All, 1, MPI_INT, world);
	  int Size_All[RankID], Size_All_Q[RankID];
	  for (int i = 0; i < RankID; i++)
	  {
		  nlocal_All_Q[i] = nlocal_All[i];
		  nlocal_All[i] = nlocal_All[i] * 3;
		  if (i == 0)
		  {
			  Size_All[i] = 0;
			  Size_All_Q[i] = 0;
		  }
		  else
		  {
			  Size_All[i] = Size_All[i - 1] + nlocal_All[i - 1];
			  Size_All_Q[i] = Size_All_Q[i - 1] + nlocal_All_Q[i - 1];
		  }
	  }

	  MPI_Request request, request_Q; MPI_Status status;

	  double X[nlocal][3], Q[nlocal];
	  for (int i = 0; i < nlocal; i++)
	  {
		  X[i][0] = x[i][0] - (boxhi[0] - Lx / 2);
		  X[i][1] = x[i][1] - (boxhi[1] - Ly / 2);
		  X[i][2] = x[i][2] - (boxhi[2] - Lz / 2);
		  //if (fabs(X[i][0]) > Lx / 2 || fabs(X[i][1]) > Ly / 2 || fabs(X[i][2]) > Lz / 2) 
			  //cout << "The exceeding point ID is " << i << endl;
		  Q[i] = q[i];
	  }

	  //cout << "   "<< X[48920][0] <<"   "<< X[48920][1] <<"    " << X[48920][2] << endl;

	  AdjustParticle_Double(X, nlocal, Lx, Ly, Lz);

	  //cout << "   " << X[48920][0] << "   " << X[48920][1] << "    " << X[48920][2] << endl;

	  double time;
	  time = MPI_Wtime();

	  MPI_Iallgatherv((double*)X, nlocal * 3, MPI_DOUBLE, (double*)AllSource, nlocal_All, Size_All, MPI_DOUBLE, world, &request);
	  MPI_Iallgatherv((double*)Q, nlocal, MPI_DOUBLE, (double*)AllQ, nlocal_All_Q, Size_All_Q, MPI_DOUBLE, world, &request_Q);

	  //Find the image charge(这里每个核只存储局部的镜像)
	  int lx = ceil((Rs - Lx / 2) / Lx), ly = ceil((Rs - Ly / 2) / Ly), lz = ceil((Rs - Lz / 2) / Lz);
	  int TotalNumber = ceil(nlocal * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);//减少存储量消耗
	  double ImageCharge[TotalNumber][4];
	  int ImageNumber;
	  SetImageCharge(ImageCharge, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), X, Q, nlocal, Rs, Lx, Ly, Lz, lx, ly, lz);

	  //Calculate near field potential which contains all the contributions in the cut-off sphere B at the monitoring points
	  double Near[Nw], Near_All[Nw];
	  //double PointSum_New[Nw][3], QuizSum_New[Nw][3];
	  //for (int i = 0; i < Nw; i++)
		  //for (int j = 0; j < 3; j++)
		  //{
			//  PointSum_New[i][j] = PointSum[i][j];
			// QuizSum_New[i][j] = QuizSum[i][j];
		  //}
	  //CalculateNearFieldAndZD(Near, ImageCharge, ImageNumber, Nw, IF_FMM_RightTerm, Lx, Ly, Lz, PointSum_New, QuizSum_New);
	  CalculateNearFieldAndZD(Near, ImageCharge, ImageNumber, Nw, IF_FMM_RightTerm, Lx, Ly, Lz, PointSum, QuizSum,tolerance);
	  MPI_Reduce(Near, Near_All, Nw, MPI_DOUBLE, MPI_SUM, 0, world);

	  //Solve The Least Square Problem
	  double C[p * p];
	  if (me == 0) { SolveLeastSquareProblem(C, A, Near_All, p, Nw); }
	  //double C[p * p], A_New[Nw][p * p - 1];
	  //for (int i = 0; i < Nw; i++)
		//  for (int j = 0; j < p * p - 1; j++)
			//  A_New[i][j] = A[i][j];
	  //if (me == 0) { SolveLeastSquareProblem(C, (double**)A_New, Near_All, p, Nw); }
	  MPI_Bcast(C, p * p, MPI_DOUBLE, 0, world);

	  MPI_Wait(&request, &status);
	  MPI_Wait(&request_Q, &status);

	  //time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for gatherv is " << time << endl;
	  //time = MPI_Wtime();

	  //Compute final force, potential and energy

	  double Energy_HSMA;
	  double Force[nlocal][3], Pot[nlocal];
	  TotalNumber = ceil(maxatom * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);//减少存储量消耗
	  double ImageCharge_All[TotalNumber][4];
	  SetImageCharge(ImageCharge_All, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), AllSource, AllQ, maxatom, Rs, Lx, Ly, Lz, lx, ly, lz);

	  //double Fibonacci_New[Np][4];
	  //for (int i = 0; i < Np; i++)
		  //for (int j = 0; j < 4; j++)
			  //Fibonacci_New[i][j] = Fibonacci[i][j];

	  //Energy_HSMA = FinalCalculateEnergyAndForce(Force, Pot, X, Q, nlocal, ImageCharge_All, ImageNumber, Fibonacci_New, (double**)QRD, (double**)QLocalRD, C, p, Fp, F, Rs, pi, IF_FMM_FinalPotential);
	  Energy_HSMA = FinalCalculateEnergyAndForce(Force, Pot, X, Q, nlocal, ImageCharge_All, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, C, p, Fp, F, Rs, pi, IF_FMM_FinalPotential,tolerance);

	  //double Energy_HSMA;
	  //double Force_All[maxatom][3], Pot_All[maxatom];
	  //Energy_HSMA = FinalCalculateEnergyAndForce(Force_All, Pot_All, AllSource, AllQ, maxatom, ImageCharge, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, C, p, Fp, F, Rs, pi, IF_FMM_FinalPotential);
	  //double Force[nlocal][3],Pot[nlocal];
	  //MPI_Reduce_scatter((double *)Force_All,(double *)Force, nlocal_All, MPI_DOUBLE, MPI_SUM,world);
	  //MPI_Reduce_scatter((double*)Pot_All, (double*)Pot, nlocal_All_Q, MPI_DOUBLE, MPI_SUM, world);

	  time = MPI_Wtime() - time;
	  //if (me == 0)cout << "Time for Final is " << time << endl;

	  scale = 1.0;
	  const double qscale = qqrd2e * scale;

	  if (me == 0) {
		  //cout << "The forces of particle 0 are" << f[0][0] << "   " << f[0][1] << "   " << f[0][2] << endl;
		  //cout << "The forces of particle 1 are" << f[1][0] << "   " << f[1][1] << "   " << f[1][2] << endl;
		  //cout << "The forces of particle 2 are" << f[2][0] << "   " << f[2][1] << "   " << f[2][2] << endl;
		  //cout << "The forces of particle 3 are" << f[3][0] << "   " << f[3][1] << "   " << f[3][2] << endl;
		  //cout << "The forces of particle 4 are" << f[4][0] << "   " << f[4][1] << "   " << f[4][2] << endl;
		  //cout << "The forces of particle 5 are" << f[5][0] << "   " << f[5][1] << "   " << f[5][2] << endl;
		  //cout << "The forces of particle 6 are" << f[6][0] << "   " << f[6][1] << "   " << f[6][2] << endl;
		  //cout << "The forces of particle 7 are" << f[7][0] << "   " << f[7][1] << "   " << f[7][2] << endl;
		  //cout << C[0] << "   " << C[1] << "    " << C[2] << endl;
		  //cout << PointSum[10][0] << "   " << PointSum[10][1] << "    " << PointSum[10][2] << endl;
		  //cout << QuizSum[10][0] << "   " << QuizSum[10][1] << "    " << QuizSum[10][2] << endl;
		  //cout << Near[0] << "   " << Near[1] << "    " << Near[2] << endl;
		  //cout << Near_All[0] << "   " << Near_All[1] << "    " << Near_All[2] << endl;
		  //cout << x[0][0] << "   " << x[0][1] << "   " << x[0][2] << endl;
		  //cout << Force[0][0] << "   " << Force[0][1] << "   " << Force[0][2] << endl;
		  //cout << Force[10000][0] << "   " << Force[10000][1] << "   " << Force[10000][2] << endl;
		  //cout <<"The potentials of particle 0 1 2 are  "<< Pot[0] << "   " << Pot[1] << "   " << Pot[2] << endl;
		  //cout << Pot[10] << "   " << Pot[11] << "   " << Pot[12] << endl;
		  //cout << "The forces of particle 0 are  " << f[0][0] << "   " << f[0][1] << "   " << f[0][2] << endl;
		  //cout << "The forces of particle 48920 are  " << f[48920][0] << "   " << f[48920][1] << "   " << f[48920][2] << endl;
	  }

	  //Assemble results
	  for (int i = 0; i < nlocal; i++)
	  {
		  f[i][0] += Force[i][0] * qscale;
		  f[i][1] += Force[i][1] * qscale;
		  f[i][2] += Force[i][2] * qscale;
	  }

	  double Energy_HSMA_All;
	  MPI_Allreduce(&Energy_HSMA, &Energy_HSMA_All, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  //energy += Energy_HSMA_All * qscale;
	  if (eflag_global) {
		  energy += Energy_HSMA_All * qscale / 2;
		  //energy += Energy_HSMA * qscale; 
	  }
	  if (evflag_atom) {
		  if (eflag_atom) {
			  for (int i = 0; i < nlocal; i++) {
				  eatom[i] += Pot[i];
				  eatom[i] *= qscale;
			  }
		  }

		  if (vflag_atom)
			  for (int i = 0; i < nlocal; i++)
				  for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
	  }

	  if (me == 0) {
		  //cout << C[0] << "   " << C[1] << "    " << C[2] << endl;
		  //cout << PointSum[10][0] << "   " << PointSum[10][1] << "    " << PointSum[10][2] << endl;
		  //cout << QuizSum[10][0] << "   " << QuizSum[10][1] << "    " << QuizSum[10][2] << endl;
		  //cout << Near[0] << "   " << Near[1] << "    " << Near[2] << endl;
		  //cout << Near_All[0] << "   " << Near_All[1] << "    " << Near_All[2] << endl;
		  //cout << "The forces of particle 0 are" << f[0][0] << "   " << f[0][1] << "   " << f[0][2] << endl;
		  //cout << "The forces of particle 0 are" << Force[0][0] << "   " << Force[0][1] << "   " << Force[0][2] << endl;
		  //cout << "The Energy is " << energy << "      " << Energy_HSMA_All * 10 / maxatom << endl;
		  //cout << "The forces of particle 48920 are  " << f[48920][0] << "   " << f[48920][1] << "   " << f[48920][2] << endl;
	  }

  }
  //cout << AllSource[1][1] << "   " << AllSource[1][2] << "    " << AllSource[1][0] << endl;
  //cout << "The information of image is "<<maxatom << "   " << ImageNumber << "    " << TotalNumber << endl;


}

double HSMA::memory_usage() 
{
  double bytes = 0.0;
  bytes += maxatom * sizeof(double);
  bytes += 3*maxatom * sizeof(double);
  return bytes;
}

double fac(double t)//calculate factorial
{
	double s;
	if (abs(t - 1) < 0.001 || abs(t) < 0.001)
		s = 1.0;
	else
	{
		s = t * fac(t - 1) + 0.00;
	}
	return s;
}

void HSMA::CalculateRDMultipoleExpansion(double* Q, int p, double x, double y, double z)
{
	Q[0] = 1.0;
	Q[1] = y / 2;
	Q[2] = (-z);
	Q[3] = (-x / 2);
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = (-(x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-z) * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	double rr = sqrt(x * x + y * y + z * z);
	t = 0;
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * (n + 0.00) / rr;
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA::CalculateLocalRDMultipoleExpansion(double* Q, int p, double x, double y, double z, double Rs)
{
	Q[0] = 1.0;
	Q[1] = y / 2;
	Q[2] = (-z);
	Q[3] = (-x / 2);
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = (-(x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-z) * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	double rr = sqrt(x * x + y * y + z * z);
	t = 0;
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * (-n - 1) * pow(Rs, 2 * n + 1.0) / pow(rr, 2 * n + 2.0);
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA::CalculateMultipoleExpansion(double* Q, int p, double x, double y, double z)
{
	Q[0] = 1.0;
	Q[1] = y / 2;
	Q[2] = (-z);
	Q[3] = (-x / 2);
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = (-(x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-z) * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}


	t = 0;//normlization     Please do not normlize after every step!! That's wrong!      
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA::CalculateZDerivativeMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	/*         先计算不求导的 开始          */
	double Qold[p * p];
	Qold[0] = 1.0;
	Qold[1] = y / 2;
	Qold[2] = (-1.0) * z;
	Qold[3] = (-1.0) * x / 2.0;
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Qold[t] = (y * Qold[n * n - 1] - x * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Qold[t] = (-(x * Qold[n * n - 1] + y * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Qold[t] = (-z) * Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Qold[t] = (-((2 * n - 1.0) * z * Qold[n * n - n + m] + (x * x + y * y + z * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}


	/*         计算不求导的 结束                       */
	Q[0] = 0.00;
	Q[1] = 0.00;
	Q[2] = -1.0;
	Q[3] = 0.00;
	t = 4;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = ((-1.0) * (x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-1.0) * z * Q[n * n - n + m] - Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (2 * n - 1.0) * Qold[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2] + (2.0 * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}

}

void HSMA::CalculateXDMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	/*         先?扑悴磺蟮?的 开?          */
	double Qold[p * p];
	Qold[0] = 1.0;
	Qold[1] = y / 2;
	Qold[2] = (-1.0) * z;
	Qold[3] = (-1.0) * x / 2.0;
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Qold[t] = (y * Qold[n * n - 1] - x * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Qold[t] = (-(x * Qold[n * n - 1] + y * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Qold[t] = (-z) * Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Qold[t] = (-((2 * n - 1.0) * z * Qold[n * n - n + m] + (x * x + y * y + z * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}


	//开?求?部分
	Q[0] = 0.00;
	Q[1] = 0.00;
	Q[2] = 0.0;
	Q[3] = -1.0 / 2.0;
	t = 4;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1] - Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = ((-1.0) * (x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1] + Qold[n * n - 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-1.0) * z * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2] + (2.0 * x) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
				//mm++;
			}
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA::CalculateYDMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	/*         先?扑悴磺蟮?的 开?          */
	double Qold[p * p];
	Qold[0] = 1.0;
	Qold[1] = y / 2;
	Qold[2] = (-1.0) * z;
	Qold[3] = (-1.0) * x / 2.0;
	int t = 4;
	int m, n;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Qold[t] = (y * Qold[n * n - 1] - x * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Qold[t] = (-(x * Qold[n * n - 1] + y * Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Qold[t] = (-z) * Qold[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Qold[t] = (-((2 * n - 1.0) * z * Qold[n * n - n + m] + (x * x + y * y + z * z) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
			}
			t++;
		}
	}

	/*         ?扑悴磺蟮?的 ?崾?                      */
		//开?求?部分

	Q[0] = 0.00;
	Q[1] = 1.0 / 2.0;
	Q[2] = 0.00;
	Q[3] = 0.00;
	t = 4;
	for (int i = 2; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			if (m == -n)
			{
				Q[t] = (y * Q[n * n - 1] - x * Q[n * n - 2 * n + 1] + Qold[n * n - 1]) / (2 * n + 0.00);
			}
			else if (m == n)
			{
				Q[t] = ((-1.0) * (x * Q[n * n - 1] + y * Q[n * n - 2 * n + 1] + Qold[n * n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				Q[t] = (-1.0) * z * Q[n * n - n + m];
			}
			else if ((n - abs(m)) > 1)
			{
				Q[t] = (-((2 * n - 1.0) * z * Q[n * n - n + m] + (x * x + y * y + z * z) * Q[n * n - 3 * n + m + 2] + (2.0 * y) * Qold[n * n - 3 * n + m + 2]) / ((n - abs(m) + 0.0) * (n + abs(m) + 0.0)));
				//mm++;
			}
			t++;
		}
	}

	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1) * (i + 1))
		{
			m = t - i - i * i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m)) * fac(n + m));
			t++;
		}
	}
}

void HSMA::SetFibonacci(double Fibonacci[][4], double F, double Fp, int Np, double Rs, double PI)
{
	double zz, zt, zs, phia, phib;
	double deltaz = 2.0 / F;
	for (int i = 0; i <= int(F); i++)
	{
		zz = -1.0 + i * deltaz;
		zt = zz + sin(PI * zz) / PI;
		zs = sqrt(1 - zt * zt);
		phia = PI * (i + 0.00) * Fp / F;
		phib = PI + phia;

		Fibonacci[2 * i][0] = cos(phia) * zs * (Rs + 0.0);
		Fibonacci[2 * i + 1][0] = cos(phib) * zs * (Rs + 0.00);
		Fibonacci[2 * i][1] = sin(phia) * zs * (Rs + 0.00);
		Fibonacci[2 * i + 1][1] = sin(phib) * zs * (Rs + 0.00);
		Fibonacci[2 * i][2] = zt * (Rs + 0.00);
		Fibonacci[2 * i + 1][2] = zt * (Rs + 0.00);
		Fibonacci[2 * i][3] = PI * deltaz * (1.0 + cos(PI * zz)) * Rs * Rs;
		Fibonacci[2 * i + 1][3] = PI * deltaz * (1.0 + cos(PI * zz)) * Rs * Rs;
	}
}

void HSMA::AdjustParticle_Float(float Particle[][3], int N, float Lx, float Ly, float Lz)//Image Correction 0~Lx
{
	for (int i = 0; i < N; i++)
	{
		if (Particle[i][0] > Lx)
		{
			Particle[i][0] = Particle[i][0] - ceil((Particle[i][0] - Lx) / (Lx + 0.00)) * Lx;
		}
		if (Particle[i][0] < 0.00)
		{
			//cout<<Particle[i].x<<"   ";
			Particle[i][0] = Particle[i][0] + ceil(fabs(Particle[i][0] / Lx)) * Lx;
			//cout<<Particle[i].x<<endl;
		}
		if (Particle[i][1] > Ly)
		{
			Particle[i][1] = Particle[i][1] - ceil((Particle[i][1] - Ly) / (Ly + 0.00)) * Ly;
		}
		if (Particle[i][1] < 0.00)
		{
			Particle[i][1] = Particle[i][1] + ceil(fabs(Particle[i][1] / Ly)) * Ly;
		}
		if (Particle[i][2] > Lz)
		{
			Particle[i][2] = Particle[i][2] - ceil((Particle[i][2] - Lz) / (Lz + 0.00)) * Lz;
		}
		if (Particle[i][2] < 0.00)
		{
			Particle[i][2] = Particle[i][2] + ceil(fabs(Particle[i][2] / Lz)) * Lz;
		}
	}
}

void HSMA::AdjustParticle_Double(double Particle[][3], int N, double Lx, double Ly, double Lz)//Image Correction -Lx/2~Lx/2
{
	for (int i = 0; i < N; i++)
	{
		if (Particle[i][0] > Lx/2)
		{
			Particle[i][0] = Particle[i][0] - ceil((Particle[i][0] - Lx/2.0) / (Lx + 0.00)) * Lx;
		}
		if (Particle[i][0] < -Lx / 2)
		{
			//cout<<Particle[i].x<<"   ";
			Particle[i][0] = Particle[i][0] + ceil(fabs((Particle[i][0] + Lx / 2.0)/Lx)) * Lx;
			//cout<<Particle[i].x<<endl;
		}
		if (Particle[i][1] > Ly/2)
		{
			Particle[i][1] = Particle[i][1] - ceil((Particle[i][1] - Ly/2.0) / (Ly + 0.00)) * Ly;
		}
		if (Particle[i][1] < -Ly / 2)
		{
			Particle[i][1] = Particle[i][1] + ceil(fabs( (Particle[i][1] + Ly / 2.0) / Ly)) * Ly;
		}
		if (Particle[i][2] > Lz/2)
		{
			Particle[i][2] = Particle[i][2] - ceil((Particle[i][2] - Lz/2.0) / (Lz + 0.00)) * Lz;
		}
		if (Particle[i][2] < -Lz / 2)
		{
			Particle[i][2] = Particle[i][2] + ceil(fabs((Particle[i][2] + Lz / 2.0)/Lz)) * Lz;
		}
	}
}

void HSMA::SetImageCharge(double ImageCharge[][4], int* ImageNumber, int TotalNumber, double Source[][3], double* Q, int NSource, double Rs, double Lx, double Ly, double Lz, int lx, int ly, int lz)
{
	int total = 0;
	int number = 0;

	for (int i = -lx; i <= lx; i++)
	{
		double CX, CY, CZ;
		for (int j = -ly; j <= ly; j++)
			for (int k = -lz; k <= lz; k++)
			{

				for (int m = 0; m < NSource; m++)
				{

					CX = Source[m][0] + i * Lx;
					CY = Source[m][1] + j * Ly;
					CZ = Source[m][2] + k * Lz;
					if (CX * CX + CY * CY + CZ * CZ <= Rs * Rs)
					{
						//#pragma omp atomic						   
						ImageCharge[number][0] = CX;
						ImageCharge[number][1] = CY;
						ImageCharge[number][2] = CZ;
						ImageCharge[number][3] = Q[m];
						number++;
					}
				}

			}
	}
	*ImageNumber = number;
}

void HSMA::CalculateNearFieldAndZD(double* Near, double ImageCharge[][4], int ImageNumber, int Nw, int IF_FMM_RightTerm, double Lx, double Ly, double Lz, double **PointSum, double **QuizSum, double tolerance)
{
	if (IF_FMM_RightTerm)//Using FMM to calculate pairwise sum
	{
		double eps = tolerance;

		/*            开始 FMM 参数设置           */
		int ns = ImageNumber;
		int nt = Nw * 2;
		double* source = (double*)malloc(3 * ns * sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double* charge = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			source[3 * i] = ImageCharge[i][0];
			source[3 * i + 1] = ImageCharge[i][1];
			source[3 * i + 2] = ImageCharge[i][2];
			charge[i] = ImageCharge[i][3];
		}

		for (int i = 0; i < Nw; i++)
		{
			target[3 * i] = PointSum[i][0];
			target[3 * i + 1] = PointSum[i][1];
			target[3 * i + 2] = PointSum[i][2];

			target[3 * Nw + 3 * i] = QuizSum[i][0];
			target[3 * Nw + 3 * i + 1] = QuizSum[i][1];
			target[3 * Nw + 3 * i + 2] = QuizSum[i][2];
		}

		double* pottarg = (double*)malloc(nt * sizeof(double));
		double* gradtarg = (double*)malloc(3 * nt * sizeof(double));

		int ier = 0;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg, &ier);		

		for (int i = 0; i < Nw; i++)
		{
			Near[i] = pottarg[Nw + i] - pottarg[i];
		}

		free(source); free(target); free(charge); free(pottarg);
		source = NULL; target = NULL; charge = NULL; pottarg = NULL;
	}
	else
	{

		double Paramet1[ImageNumber];
		for (int i = 0; i < ImageNumber; i++)
		{
			Paramet1[i] = ImageCharge[i][3];
		}

		double Paramet[ImageNumber];
		memcpy(Paramet, Paramet1, sizeof(double) * ImageNumber);


		double Image[ImageNumber][4];
		memcpy((double*)Image, (double*)ImageCharge, sizeof(double) * 4 * ImageNumber);


	#pragma omp parallel shared(Near,Paramet,Image) 
		{

			double pottarg, fldtarg, pottarg2, fldtarg2;
			double deltax, deltay, delta;
			double deltaz;

			double deltax1, deltay1, deltaz1;

			double delta1;
		#pragma omp for 
			for (int i0 = 0; i0 < Nw; i0 += 20)
			{
				//#pragma vector nontemporal
				for (int i = i0; i < (i0 + 20 < Nw ? i0 + 20 : Nw); i++)
				{
					pottarg = 0.00;
					fldtarg = 0.00;
					pottarg2 = 0.00;
					fldtarg2 = 0.00;
					double Copy1 = PointSum[i][0], Copy2 = PointSum[i][1], Copy3 = PointSum[i][2], Copy4 = QuizSum[i][0], Copy5 = QuizSum[i][1], Copy6 = QuizSum[i][2];
					for (int j = 0; j < ImageNumber; j++)
					{
						deltax = Copy1 - Image[j][0];
						deltay = Copy2 - Image[j][1];
						deltaz = Copy3 - Image[j][2];

						double Para = Paramet[j];

						delta = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
						pottarg = pottarg + Para / delta;

						//deltaz1=-lz-Image[j][2]; 
						deltax1 = Copy4 - Image[j][0];
						deltay1 = Copy5 - Image[j][1];
						deltaz1 = Copy6 - Image[j][2];

						delta1 = sqrt(deltax1 * deltax1 + deltay1 * deltay1 + deltaz1 * deltaz1);
						pottarg2 = pottarg2 + Para / delta1;
					}

					Near[i] = pottarg2 - pottarg;
				}
			}


		}

	}

}

void HSMA::SolveLeastSquareProblem(double* C, double** A, double* Near, int p, int Nw)
{

	int rowAT, columnATA, columnAT;
	double alpha, beta;
	rowAT = p * p - 1; columnAT = Nw; columnATA = p * p - 1;
	double MatrixAT[rowAT * columnAT], MatrixATA[rowAT * columnATA], BB[Nw], ATB[p * p - 1], INV_ATA_ATB[p * p - 1];

	alpha = 1.0; beta = 0.00;
	

	for (int i = 0; i < rowAT; i++)
		for (int j = 0; j < columnAT; j++)
		{
			MatrixAT[i * columnAT + j] = A[j][i];
		}

	for (int i = 0; i < (rowAT * columnATA); i++) {
		MatrixATA[i] = 0.0;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowAT, columnATA, columnAT, alpha, MatrixAT, columnAT, MatrixAT, columnAT, beta, MatrixATA, columnATA);
	
	int InfoHelp;
	int VectorHelp[columnATA];
	//int * VectorHelp;
	//VectorHelp = (int*)mkl_malloc(columnATA * sizeof(int), 64);
	
	for (int i = 0; i < columnATA; i++)
		VectorHelp[i] = 0;
	InfoHelp = LAPACKE_dgetrf(CblasRowMajor, columnATA, columnATA, MatrixATA, columnATA, VectorHelp);
	InfoHelp = LAPACKE_dgetri(CblasRowMajor, columnATA, MatrixATA, columnATA, VectorHelp);

	for (int i = 0; i < columnAT; i++)
	{
		BB[i] = Near[i];
	}
	for (int i = 0; i < columnATA; i++)
	{
		ATB[i] = 0.00;
		INV_ATA_ATB[i] = 0.00;
	}

	cblas_dgemv(CblasRowMajor, CblasNoTrans, rowAT, columnAT, alpha, MatrixAT, columnAT, BB, 1, beta, ATB, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, columnATA, columnATA, alpha, MatrixATA, columnATA, ATB, 1, beta, INV_ATA_ATB, 1);


	C[0] = 0.00;
	for (int i = 0; i < rowAT; i++)
	{
		C[i + 1] = INV_ATA_ATB[i];
	}
}

double HSMA::FinalCalculateEnergyAndForce(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][4], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance)
{
	if (!IF_FMM_FinalPotential)
	{

		double EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
		for (int i = 0; i < NSource; i++)
		{
			EF[i] = 0.00; EFX[i] = 0.00; EFY[i] = 0.00; EFZ[i] = 0.00;
		}
		//double startTime = omp_get_wtime();
	#pragma omp parallel shared(Source,EF,EFX,EFY,EFZ,p,C)   
		{
			double QF[p * p], QFX[p * p], QFY[p * p], QFZ[p * p];
			double CC[p * p];
			memcpy(CC, C, sizeof(double) * p * p);

	        #pragma omp for //reduction(+:EF,EFX,EFY,EFZ)
			for (int i = 0; i < NSource; i++)
			{
				CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
				CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
				CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
				CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
				//EF[i]=0.00;EFX[i]=0.00;EFY[i]=0.00;EFZ[i]=0.00;
				double a = 0.00, b = 0.00, c = 0.00, d = 0.00;
				for (int j = 0; j < p * p; j++)
				{
					a = a + QF[j] * CC[j];
					b = b + QFX[j] * CC[j];
					c = c + QFY[j] * CC[j];
					d = d + QFZ[j] * CC[j];
				}
				EF[i] = EF[i] + a;
				EFX[i] = EFX[i] + b;
				EFY[i] = EFY[i] + c;
				EFZ[i] = EFZ[i] + d;
			}
		}
		//double endTime=omp_get_wtime();
		//cout<<"The first part of Update Step II Cost Time:  "<<endTime-startTime<<endl<<endl;
		
		double startTime = omp_get_wtime();
		double EN[NSource],ENX[NSource],ENY[NSource],ENZ[NSource];
		for(int i=0;i<NSource;i++)
		{
			EN[i]=0.00; ENX[i]=0.00; ENY[i]=0.00; ENZ[i]=0.00;
		}

		double QImage[ImageNumber];
		for(int j=0;j<ImageNumber;j++)
		{
			QImage[j]=ImageCharge[j][3];
		}

		double Image[ImageNumber][4];
		memcpy((double *)Image,(double *)ImageCharge,sizeof(double)*4*ImageNumber);

		#pragma omp parallel
		{
			double deltax,deltay,deltaz,delta;

			#pragma omp for
			for(int i=0;i<NSource;i++)
			{
				double a=0.00,b=0.00,c=0.00,d=0.00;
				for(int j=0;j<ImageNumber;j++)
				{
					deltax=(Image[j][0]-Source[i][0]);
					deltay=(Image[j][1]-Source[i][1]);
					deltaz=(Image[j][2]-Source[i][2]);

				delta=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
				if(!((fabs(deltax)<1.0e-13)&&(fabs(deltay)<1.0e-13)&&(fabs(deltaz)<1.0e-13))  )
					{
						a=a+QImage[j]/delta;
						b=b+QImage[j]*(deltax)/(delta*delta*delta);
						c=c+QImage[j]*(deltay)/(delta*delta*delta);
						d=d+QImage[j]*(deltaz)/(delta*delta*delta);
					}
				}
				EN[i]=a;
				ENX[i]=b;
				ENY[i]=c;
				ENZ[i]=d;
			}
		}

		//double endTime=omp_get_wtime();
		//cout<<"The first part of Update Step II Cost Time:  "<<endTime-startTime<<endl<<endl;
		double Energy = 0.00;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = EN[i] + EF[i];
			Energy = Energy + Q[i] * Pot[i];
			Force[i][0] = -(EFX[i] + ENX[i]) * Q[i];
			Force[i][1] = -(EFY[i] + ENY[i]) * Q[i];
			Force[i][2] = -(EFZ[i] + ENZ[i]) * Q[i];

		}
		if (me == 0) {
			//cout << "开始   " <<setprecision(16)<< EN[1]<<"   "<<EF[1]<<"    "<<10*(EN[1] + EF[1]) << endl;
			//cout << ENX[1] + EFX[1] << endl;
			//cout << ENY[1] + EFY[1] << endl;
			//cout << ENZ[1] + EFZ[1] << endl;
		}

		return Energy;
	}
	else
	{
		double eps = tolerance;

		int ns = ImageNumber;
		int nt = NSource;

		double *source = (double *)malloc(3*ns*sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double *charge = (double *)malloc(ns*sizeof(double));

		for(int i=0;i<ns;i++)
		{
			source[3*i]=ImageCharge[i][0];
			source[3*i+1]=ImageCharge[i][1];
			source[3*i+2]=ImageCharge[i][2];
			charge[i]=ImageCharge[i][3];	
		}

		for (int i = 0; i < nt; i++)
		{
			target[3 * i] = Source[i][0];
			target[3 * i + 1] = Source[i][1];
			target[3 * i + 2] = Source[i][2];
		}

		double *pottarg = (double *)malloc(nt*sizeof(double));
		double *gradtarg = (double *)malloc(3*nt*sizeof(double));

		int ier1;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg,&ier1);

		int KL = int(2 * F + 2);
		/*             BEGIN HSMA ALGORITHM                          */
		double EF[KL];
		double CenterPara;

		//#pragma omp parallel for
		for (int i = 0; i < KL; i++)
		{
			CenterPara = 0.00;
			for (int j = 0; j < p * p; j++)
			{
				CenterPara = CenterPara + (1 / (4 * PI)) * (QRD[i][j] - QLocalRD[i][j]) * C[j];
			}
			EF[i] = CenterPara * Fibonacci[i][3];
		}

		/*            开始 FMM 参数设置           */
		ns = 2 * F + 2;
		double* sourceF = (double*)malloc(3 * ns * sizeof(double));
		double* chargeF = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			sourceF[3 * i] = Fibonacci[i][0];
			sourceF[3 * i + 1] = Fibonacci[i][1];
			sourceF[3 * i + 2] = Fibonacci[i][2];
			chargeF[i] = EF[i];
		}

		double* pottargF = (double*)malloc(nt * sizeof(double));
		double* gradtargF = (double*)malloc(3 * nt * sizeof(double));

		int ier;
		lfmm3d_t_c_g_(&eps, &ns, sourceF, chargeF, &nt, target, pottargF, gradtargF,&ier);

		/*	Final Summation	*/
		double Energy = 0.00;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = pottargF[i]+ pottarg[i];
			Energy = Energy + Pot[i] * Q[i];
			Force[i][0] = - (gradtargF[3 * i]+ gradtarg[3 * i]) * Q[i];
			Force[i][1] = - (gradtargF[3 * i + 1]+ gradtarg[3 * i+1]) * Q[i];
			Force[i][2] = - (gradtargF[3 * i + 2]+ gradtarg[3 * i+2]) * Q[i];
		}

		free(target);
		target = NULL;
		free(sourceF); free(chargeF); free(pottargF); free(gradtargF);
		sourceF = NULL; chargeF = NULL; pottargF = NULL; gradtargF = NULL;
		 
		return Energy;
	}

}

void HSMA::CalculateNearFieldAndZD_Single(double* Near, double ImageCharge[][4], int ImageNumber, int Nw, int IF_FMM_RightTerm, double Lx, double Ly, double Lz, double **PointSum, double **QuizSum, double Source[][3], double Force[][3], double* Pot, int NSource, double* Q, double tolerance)
{
	if (IF_FMM_RightTerm)//Using FMM to calculate pairwise sum
	{
		double eps = tolerance;

		/*            开始 FMM 参数设置           */
		int ns = ImageNumber;
		int nt = Nw * 2 + NSource;
		double* source = (double*)malloc(3 * ns * sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double* charge = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			source[3 * i] = ImageCharge[i][0];
			source[3 * i + 1] = ImageCharge[i][1];
			source[3 * i + 2] = ImageCharge[i][2];
			charge[i] = ImageCharge[i][3];
		}

		for (int i = 0; i < Nw; i++)
		{
			target[3 * i] = PointSum[i][0];
			target[3 * i + 1] = PointSum[i][1];
			target[3 * i + 2] = PointSum[i][2];

			target[3 * Nw + 3 * i] = QuizSum[i][0];
			target[3 * Nw + 3 * i + 1] = QuizSum[i][1];
			target[3 * Nw + 3 * i + 2] = QuizSum[i][2];
		}

		for (int i = 0; i < NSource; i++)
		{
			target[6 * Nw + 3 * i] = Source[i][0];
			target[6 * Nw + 3 * i + 1] = Source[i][1];
			target[6 * Nw + 3 * i + 2] = Source[i][2];
		}

		double* pottarg = (double*)malloc(nt * sizeof(double));
		double* gradtarg = (double*)malloc(3 * nt * sizeof(double));

		//double startTime = omp_get_wtime();
		int ier;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg,&ier);
		//double endTime = omp_get_wtime();
		//cout<<endl<<"Evaluate Cost Time="<<endTime-startTime<<endl<<endl<<endl;		

		for (int i = 0; i < Nw; i++)
		{
			Near[i] = pottarg[Nw + i] - pottarg[i];
		}
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = pottarg[2 * Nw + i];
			Force[i][0] = -gradtarg[6 * Nw + 3 * i] * Q[i];
			Force[i][1] = -gradtarg[6 * Nw + 3 * i + 1] * Q[i];
			Force[i][2] = -gradtarg[6 * Nw + 3 * i + 2] * Q[i];
		}


		free(source); free(target); free(charge); free(pottarg);
		source = NULL; target = NULL; charge = NULL; pottarg = NULL;
	}
	else
	{

		double Paramet1[ImageNumber];
		for (int i = 0; i < ImageNumber; i++)
		{
			Paramet1[i] = ImageCharge[i][3];
		}

		double Paramet[ImageNumber];
		memcpy(Paramet, Paramet1, sizeof(double) * ImageNumber);


		double Image[ImageNumber][4];
		memcpy((double*)Image, (double*)ImageCharge, sizeof(double) * 4 * ImageNumber);


	#pragma omp parallel shared(Near,Paramet,Image,Source,Force,Pot,Q,NSource) 
		{

			double pottarg, fldtarg, pottarg2, fldtarg2;
			double deltax, deltay, delta;
			double deltaz;

			double lz = Lz / 2.0;

			double deltax1, deltay1, deltaz1;


			double delta1;
	#pragma omp for 
			for (int i0 = 0; i0 < Nw; i0 += 20)
			{
				//#pragma vector nontemporal
				for (int i = i0; i < (i0 + 20 < Nw ? i0 + 20 : Nw); i++)
				{
					pottarg = 0.00;
					fldtarg = 0.00;
					pottarg2 = 0.00;
					fldtarg2 = 0.00;
					double Copy1 = PointSum[i][0], Copy2 = PointSum[i][1], Copy3 = PointSum[i][2], Copy4 = QuizSum[i][0], Copy5 = QuizSum[i][1], Copy6 = QuizSum[i][2];
					for (int j = 0; j < ImageNumber; j++)
					{
						deltax = Copy1 - Image[j][0];
						deltay = Copy2 - Image[j][1];
						deltaz = Copy3 - Image[j][2];

						double Para = Paramet[j];

						delta = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
						pottarg = pottarg + Para / delta;

						//deltaz1=-lz-Image[j][2]; 
						deltax1 = Copy4 - Image[j][0];
						deltay1 = Copy5 - Image[j][1];
						deltaz1 = Copy6 - Image[j][2];

						delta1 = sqrt(deltax1 * deltax1 + deltay1 * deltay1 + deltaz1 * deltaz1);
						pottarg2 = pottarg2 + Para / delta1;
					}


					Near[i] = pottarg2 - pottarg;
				}
			}
	#pragma omp for 
			for (int i0 = 0; i0 < NSource; i0 += 20)
			{
				for (int i = i0; i < (i0 + 20 < NSource ? i0 + 20 : NSource); i++)
				{
					double a = 0.00, b = 0.00, c = 0.00, d = 0.00;
					for (int j = 0; j < ImageNumber; j++)
					{
						deltax = (Image[j][0] - Source[i][0]);
						deltay = (Image[j][1] - Source[i][1]);
						deltaz = (Image[j][2] - Source[i][2]);
						double QImage = Image[j][3];
						delta = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
						if (!((fabs(deltax) < 1.0e-13) && (fabs(deltay) < 1.0e-13) && (fabs(deltaz) < 1.0e-13)))
						{
							a = a + QImage / delta;
							b = b + QImage * (deltax) / (delta * delta * delta);
							c = c + QImage * (deltay) / (delta * delta * delta);
							d = d + QImage * (deltaz) / (delta * delta * delta);
						}
					}
					Pot[i] = a;
					Force[i][0] = -b * Q[i];
					Force[i][1] = -c * Q[i];
					Force[i][2] = -d * Q[i];
				}
			}


		}

	}

}

double HSMA::FinalCalculateEnergyAndForce_Single(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][4], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance)
{
	if (!IF_FMM_FinalPotential)
	{

		double EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
		for (int i = 0; i < NSource; i++)
		{
			EF[i] = 0.00; EFX[i] = 0.00; EFY[i] = 0.00; EFZ[i] = 0.00;
		}
		//double startTime = omp_get_wtime();
	#pragma omp parallel shared(Source,EF,EFX,EFY,EFZ,p,C)   
		{
			double QF[p * p], QFX[p * p], QFY[p * p], QFZ[p * p];
			double CC[p * p];
			memcpy(CC, C, sizeof(double) * p * p);

	#pragma omp for //reduction(+:EF,EFX,EFY,EFZ)
			for (int i = 0; i < NSource; i++)
			{
				CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
				CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
				CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
				CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
				//EF[i]=0.00;EFX[i]=0.00;EFY[i]=0.00;EFZ[i]=0.00;
				double a = 0.00, b = 0.00, c = 0.00, d = 0.00;
				for (int j = 0; j < p * p; j++)
				{
					a = a + QF[j] * CC[j];
					b = b + QFX[j] * CC[j];
					c = c + QFY[j] * CC[j];
					d = d + QFZ[j] * CC[j];
				}
				EF[i] = EF[i] + a;
				EFX[i] = EFX[i] + b;
				EFY[i] = EFY[i] + c;
				EFZ[i] = EFZ[i] + d;
			}
		}

		double Energy = 0.00;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = Pot[i] + EF[i];
			Energy = Energy + Q[i] * Pot[i];
			Force[i][0] = Force[i][0] - (EFX[i]) * Q[i];
			Force[i][1] = Force[i][1] - (EFY[i]) * Q[i];
			Force[i][2] = Force[i][2] - (EFZ[i]) * Q[i];
		}
		return Energy;
	}
	else
	{

		double eps = tolerance;

		int ns = ImageNumber;
		int nt = NSource;

		double* target = (double*)malloc(3 * nt * sizeof(double));

		for (int i = 0; i < nt; i++)
		{
			target[3 * i] = Source[i][0];
			target[3 * i + 1] = Source[i][1];
			target[3 * i + 2] = Source[i][2];
		}

		int KL = int(2 * F + 2);
		/*             BEGIN HSMA ALGORITHM                          */
		double EF[KL];
		double CenterPara;

		//#pragma omp parallel for
		for (int i = 0; i < KL; i++)
		{
			CenterPara = 0.00;
			for (int j = 0; j < p * p; j++)
			{
				CenterPara = CenterPara + (1 / (4 * PI)) * (QRD[i][j] - QLocalRD[i][j]) * C[j];//*((double*)QRD + i * p * p + j) - *((double*)QLocalRD + i * p * p + j)
			}
			EF[i] = CenterPara * Fibonacci[i][3];
		}

		/*            开始 FMM 参数设置           */
		ns = 2 * F + 2;
		double* sourceF = (double*)malloc(3 * ns * sizeof(double));
		double* chargeF = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			sourceF[3 * i] = Fibonacci[i][0];
			sourceF[3 * i + 1] = Fibonacci[i][1];
			sourceF[3 * i + 2] = Fibonacci[i][2];
			chargeF[i] = EF[i];
		}

		double* pottargF = (double*)malloc(nt * sizeof(double));
		double* gradtargF = (double*)malloc(3 * nt * sizeof(double));

		int ier;
		lfmm3d_t_c_g_(&eps, &ns, sourceF, chargeF, &nt, target, pottargF, gradtargF,&ier);

		/*	Final Summation	*/
		double Energy = 0.00;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = Pot[i] + pottargF[i];
			Energy = Energy + Pot[i] * Q[i];
			Force[i][0] = Force[i][0] - (gradtargF[3 * i]) * Q[i];
			Force[i][1] = Force[i][1] - (gradtargF[3 * i + 1]) * Q[i];
			Force[i][2] = Force[i][2] - (gradtargF[3 * i + 2]) * Q[i];
		}

		free(target);
		target = NULL;
		free(sourceF); free(chargeF); free(pottargF); free(gradtargF);
		sourceF = NULL; chargeF = NULL; pottargF = NULL; gradtargF = NULL;

		return Energy;
	}

}

