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

/* ----------------------------------------------------------------------
   Contributing author: Jiuyang Liang (liangjiuyang@sjtu.edu.cn)
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
#include "update.h"
#include "hsma2d.h"
#include<iostream>
#include <sstream>
#include <fstream>
#include"mkl.h"
#include<omp.h>
#include<iomanip>
#include <immintrin.h> 
#include "complex.h"

extern "C" {void lfmm3d_t_c_g_(double* eps, int* nsource, double* source, double* charge, int* nt, double* targ, double* pottarg, double* gradtarg, int* ier); }
extern int fab(int n);
extern int isfab(int m);

using namespace LAMMPS_NS;
using namespace std;

HSMA2D::HSMA2D(LAMMPS* lmp) : KSpace(lmp)
{
	maxatom = atom->natoms;
	MPI_Comm_rank(world, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &RankID);
	Lx = domain->xprd;
	Ly = domain->yprd;
	Lz = domain->zprd * slab_volfactor;
	PI = 3.141592653589793;
}

/* ---------------------------------------------------------------------- */
void HSMA2D::settings(int narg, char** arg)
{
	if (narg != 10) error->all(FLERR, "Illegal kspace_style HSMA2D command");
	tolerance = fabs(utils::numeric(FLERR, arg[0], false, lmp));

	Lambda = utils::numeric(FLERR, arg[1], false, lmp);
	Gamma= utils::numeric(FLERR, arg[2], false, lmp);
	p = utils::numeric(FLERR, arg[3], false, lmp);
	Nw = utils::numeric(FLERR, arg[4], false, lmp);
	w = utils::numeric(FLERR, arg[5], false, lmp);
	Fp = utils::numeric(FLERR, arg[6], false, lmp);
	F = utils::numeric(FLERR, arg[7], false, lmp);
	IF_FMM_RightTerm = utils::numeric(FLERR, arg[8], false, lmp);
	IF_FMM_FinalPotential = utils::numeric(FLERR, arg[9], false, lmp);

	if (Lambda < 0)
		error->all(FLERR, "Lambda shoule be >0.");
	if (Lambda > 20 || Lambda < 0.2) {
		error->warning(FLERR, fmt::format("The Lambda is too big/small! Please use an approximate range of Lambda. Set Lambda to the default value.",
			update->ntimestep));
		Lambda = 1.5;
	}
	if (p < 1)
		error->all(FLERR, "p shoule be >=1.");
	if (p > 50) {
		error->warning(FLERR, fmt::format("The p is too big! Please use a smaller p. Set p to the default value.",
			update->ntimestep));
		p = 6;
	}
	if ((IF_FMM_RightTerm) != 0 && (IF_FMM_RightTerm != 1))
		error->all(FLERR, "Using wrong value of IF_FMM_RightTerm.");
	if ((IF_FMM_FinalPotential) != 0 && (IF_FMM_FinalPotential != 1))
		error->all(FLERR, "Using wrong value of IF_FMM_FinalPotential.");
	if (Fp > F || Fp < 0 || F < 0 || isfab(int(Fp)) != 1 || isfab(int(F)) != 1) {
		error->warning(FLERR, fmt::format("Using wrong value of Fp and F. Set Fp and F to the default value.",
			update->ntimestep));
		Fp = 89.0;
		F = 144.0;
	}
	if (Nw <= 5||Nw>=150)
	{
		error->warning(FLERR, fmt::format("Nw is too small/big! Set Nw to the default value.",
			update->ntimestep));
		Nw = 40.0;
	}
	if (w <= 0)
	{
		error->all(FLERR, "Using wrong value of w.");
	}
	if (domain->dimension == 3)
		error->all(FLERR, "Cannot use HSMA2D in a 3d simulation. Please use HSMA3D instead.");
}

/* ---------------------------------------------------------------------- */
void HSMA2D::init()
{
	printf("Setting up HSMA2D implemented by Jiuyang Liang (Release 1.0.0)\n");

	Step = 0;
	Time = new float[2000];

	Gauss = new double[Nw];
	Weight = new double[Nw];
	Gauss_Grule3D(Gauss, Weight, Nw);
	S = ceil(sqrt(w));
	double ar[S][Nw];
	for (int tt = 0; tt < S; tt++)
	{
		for (int i = 0; i < Nw; i++)
		{
			ar[tt][i] = (((Gauss[i] + 1) / 2.0) * 2.0 - S + (tt) * 2.0) / S;
		}
	}

	IntegralTop = new double** [S * Nw];
	IntegralDown = new double** [S * Nw];
	for (int i = 0; i < S * Nw; i++)
	{
		IntegralTop[i] = new double* [S * Nw];
		IntegralDown[i] = new double* [S * Nw];
	}
	for (int i = 0; i < S * Nw; i++)
	{
		for (int j = 0; j < S * Nw; j++)
		{
			IntegralTop[i][j] = new double [4];
			IntegralDown[i][j] = new double [4];
		}
	}
	for (int m = 0; m < S; m++)
		for (int i = 0; i < Nw; i++)
			for (int n = 0; n < S; n++)
				for (int j = 0; j < Nw; j++)
				{
					IntegralTop[m * Nw + i][n * Nw + j][0] = ar[m][i] * Lx / 2.0;
					IntegralTop[m * Nw + i][n * Nw + j][1] = ar[n][j] * Ly / 2.0;
					IntegralTop[m * Nw + i][n * Nw + j][2] = Lz / 2;
					IntegralTop[m * Nw + i][n * Nw + j][3] = Weight[i] * Weight[j] * Lx * Ly / (4.0 * S * S);
					IntegralDown[m * Nw + i][n * Nw + j][0] = ar[m][i] * Lx / 2.0;
					IntegralDown[m * Nw + i][n * Nw + j][1] = ar[n][j] * Ly / 2.0;
					IntegralDown[m * Nw + i][n * Nw + j][2] = -Lz / 2;
					IntegralDown[m * Nw + i][n * Nw + j][3] = Weight[i] * Weight[j] * Lx * Ly / (4.0 * S * S);
				}
	 
	Rs = Lambda * sqrt((Lx / 2) * (Lx / 2) + (Ly / 2) * (Ly / 2) + (Lz / 2) * (Lz / 2));
	PI = 3.141592653589793;
	mu = 2.0 * PI / Lx;
	EU = (1 - Gamma) / (1 + Gamma);

	IndexReal = floor((p * p - 1) / 2.0);
	IndexImag = ceil((p * p - 1) / 2.0); 
	NJKBound = floor(2.0 * p / 3.0);
	NJK = (2 * NJKBound + 1) * (2 * NJKBound + 1);
	LeftTermReal = new double* [2 * NJK];
	LeftTermImag = new double* [2 * NJK];
	for (int i = 0; i < 2 * NJK; i++)
	{
		LeftTermReal[i] = new double[IndexReal];
		LeftTermImag[i] = new double[IndexImag];
	}
	for (int i = 0; i < 2 * NJK; i++)
	{
		for (int j = 0; j < IndexReal; j++)
		{
			LeftTermReal[i][j] = 0.00;
		}
		for (int j = 0; j < IndexImag; j++)
		{
			LeftTermImag[i][j] = 0.00;
		}
	}
	ConstructLeftTerm((double**)LeftTermReal, (double**)LeftTermImag, (double***)IntegralTop, (double***)IntegralDown, Nw, NJKBound, p, EU, mu, Lx, Ly, Lz, S);//Construct the left term.

	 //Construct Fibonacci Points And There Local/Multipole Expansion.  Fibonacci[Np][4] storage the positions and the weights.
	Np = int(2 * F + 2.0);
	Fibonacci = new double* [Np];
	for (int i = 0; i < Np; i++)
	{
		Fibonacci[i] = new double[4];
	}
	QRD = new double* [Np]; QLocalRD = new double* [Np];
	for (int i = 0; i < Np; i++)
	{
		QRD[i] = new double[p * p];
		QLocalRD[i]= new double[p * p];
	}
	SetFibonacci(Fibonacci, F, Fp, Np, Rs, PI);
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

	TopNear = new double* [S * Nw];
	TopZDNear = new double* [S * Nw];
	DownNear = new double* [S * Nw];
	DownZDNear = new double* [S * Nw];
	for (int i = 0; i < S * Nw; i++)
	{
		TopNear[i] = new double[S * Nw];
		TopZDNear[i] = new double[S * Nw];
		DownNear[i] = new double[S * Nw];
		DownZDNear[i] = new double[S * Nw];
	}

	AR = new double[S * Nw];
	for (int i = 0; i < S; i++)
	{
		for (int j = 0; j < Nw; j++)
		{
			AR[i * Nw + j] = ar[i][j];
		}
	}

	cout << Lx << "=Lx   " << Ly << "=Ly    " << Lz << "=Lz   " << Lambda << "=Lambda    " <<Gamma<<"=Gamma    "<<w<<"=w    "<< p << "=p   " << Nw << "=Nw    " << Fp << "=Fp   " << F << "=F   " << IF_FMM_RightTerm << "=IF_FMM_RightTerm   " << IF_FMM_FinalPotential << "=IF_FMM_FinalPotential   " << endl;
}

void HSMA2D::compute(int eflag, int vflag)
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

	double time;
	time = MPI_Wtime();

	//Set interfaces
	double** x = atom->x;
	double* q = atom->q;
	int nlocal = atom->nlocal;
	qqrd2e = force->qqrd2e;
	double** f = atom->f;
	double boxlo[3] = { domain->boxlo[0] ,domain->boxlo[1],domain->boxlo[2] };
	double boxhi[3] = { domain->boxhi[0] ,domain->boxhi[1],domain->boxhi[2] };

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

		int lx = ceil((Rs - Lx / 2) / Lx), ly = ceil((Rs - Ly / 2) / Ly), lz = ceil((Rs - Lz / 2) / Lz);
		int TotalNumber = nlocal * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1);
		double ImageCharge[TotalNumber][5];
		int ImageNumber;
		SetImageCharge(ImageCharge, &ImageNumber, (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1), X, Q, nlocal, Rs, Lx, Ly, Lz, lx, ly, lz);

		/*Calculate Near Field and ZD*/
		CalculateNearFieldAndZD((double**)TopNear, (double**)TopZDNear, (double**)DownNear, (double**)DownZDNear, ImageCharge, ImageNumber, (double***)IntegralTop, (double***)IntegralDown, Nw, Gamma, IF_FMM_RightTerm, S, AR, Lx, Ly, Lz, tolerance);

		/*Construct the right term by using the DTN condition*/
		double RightTermReal[2 * NJK], RightTermImag[2 * NJK];
		ConstructRightTerm(RightTermReal, RightTermImag, (double**)TopNear, (double**)TopZDNear, (double**)DownNear, (double**)DownZDNear, (double***)IntegralTop, (double***)IntegralDown, Nw, NJKBound, EU, mu, Lx, Ly, Lz, S, tolerance);

		//Solve The Least Square Problem
		double C[p * p];
		SolveLeastSquareProblem(C, (double**)LeftTermReal, (double**)LeftTermImag, RightTermReal, RightTermImag, p, 2 * NJK);

		/*Calculate the potential, energyand force of the source charges*/
		double Energy_HSMA;
		Energy_HSMA = FinalCalculateEnergyAndForce(Force, Pot, X, Q, nlocal, ImageCharge, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, Gamma, C, p, Fp, F, Rs, PI, IF_FMM_FinalPotential, tolerance);

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
	}
	else if (RankID > 1) {
		double AllSource[maxatom][3], AllQ[maxatom];
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
			Q[i] = q[i];
		}
		AdjustParticle_Double(X, nlocal, Lx, Ly, Lz);
		MPI_Iallgatherv((double*)X, nlocal * 3, MPI_DOUBLE, (double*)AllSource, nlocal_All, Size_All, MPI_DOUBLE, world, &request);
		MPI_Iallgatherv((double*)Q, nlocal, MPI_DOUBLE, (double*)AllQ, nlocal_All_Q, Size_All_Q, MPI_DOUBLE, world, &request_Q);
        
		int lx = ceil((Rs - Lx / 2.0) / Lx), ly = ceil((Rs - Ly / 2.0) / Ly), lz = ceil((Rs - Lz / 2.0) / Lz);
		int TotalNumber = nlocal * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1);
		double ImageCharge[TotalNumber][5];
		int ImageNumber;
		SetImageCharge(ImageCharge, &ImageNumber, (2 * lx + 1)* (2 * ly + 1)* (2 * lz + 1), X, Q, nlocal, Rs, Lx, Ly, Lz, lx, ly, lz);

		/*Calculate Near Field and ZD*/
		CalculateNearFieldAndZD((double**)TopNear, (double**)TopZDNear, (double**)DownNear, (double**)DownZDNear, ImageCharge, ImageNumber, (double***)IntegralTop, (double***)IntegralDown, Nw, Gamma, IF_FMM_RightTerm, S, AR, Lx, Ly, Lz, tolerance);
		
		/*Construct the right term by using the DTN condition*/
		double RightTermReal[2 * NJK], RightTermImag[2 * NJK];
		ConstructRightTerm(RightTermReal, RightTermImag, (double**)TopNear, (double**)TopZDNear, (double**)DownNear, (double**)DownZDNear, (double***)IntegralTop, (double***)IntegralDown, Nw, NJKBound, EU, mu, Lx, Ly, Lz, S, tolerance);
		double RightTermReal_All[2 * NJK], RightTermImag_All[2 * NJK];
		MPI_Reduce(RightTermReal, RightTermReal_All, 2 * NJK, MPI_DOUBLE, MPI_SUM, 0, world);
		MPI_Reduce(RightTermImag, RightTermImag_All, 2 * NJK, MPI_DOUBLE, MPI_SUM, 0, world);

		//Solve The Least Square Problem
		double C[p * p];
		if (me == 0) {
			SolveLeastSquareProblem(C, (double**)LeftTermReal, (double**)LeftTermImag, RightTermReal, RightTermImag, p, 2 * NJK);
	    }
		MPI_Bcast(C, p* p, MPI_DOUBLE, 0, world);
		MPI_Wait(&request, &status);
		MPI_Wait(&request_Q, &status);

		double Energy_HSMA;
		double Force[nlocal][3], Pot[nlocal];
		TotalNumber = ceil(maxatom * (2 * lx + 1) * (2 * ly + 1) * (2 * lz + 1) / 2);//���ٴ洢������
		double ImageCharge_All[TotalNumber][5];
		SetImageCharge(ImageCharge_All, &ImageNumber, (2 * lx + 1)* (2 * ly + 1)* (2 * lz + 1), AllSource, AllQ, maxatom, Rs, Lx, Ly, Lz, lx, ly, lz);
        
		/*Calculate the potential, energyand force of the source charges*/
		Energy_HSMA = FinalCalculateEnergyAndForce(Force, Pot, X, Q, nlocal, ImageCharge_All, ImageNumber, Fibonacci, (double**)QRD, (double**)QLocalRD, Gamma,C, p, Fp, F, Rs, PI, IF_FMM_FinalPotential, tolerance);
        
		scale = 1.0;
		const double qscale = qqrd2e * scale;

		//Assemble results
		for (int i = 0; i < nlocal; i++)
		{
			f[i][0] += Force[i][0] * qscale;
			f[i][1] += Force[i][1] * qscale;
			f[i][2] += Force[i][2] * qscale;
		}

		double Energy_HSMA_All;
		MPI_Allreduce(&Energy_HSMA, &Energy_HSMA_All, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if (eflag_global) {
			energy += Energy_HSMA_All * qscale / 2;
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
	}

	time = MPI_Wtime() - time;
	int TotalStep = 100;
	if (Step < TotalStep) { Time[Step] = time * 1000; }
	Step++;
	if (me == 0 && Step == TotalStep) {
		ofstream outfile;
		outfile.open("./Time_HSMA2D.txt");
		for (int i = 0; i < TotalStep; i++)
			outfile << Time[i] << endl;
		outfile.close();
	}
}

double HSMA2D::memory_usage()
{
	double bytes = 0.0;
	bytes += maxatom * sizeof(double);
	bytes += 3 * maxatom * sizeof(double);
	return bytes;
}

void HSMA2D::Gauss_Grule3D(double* Point, double* Weight, int N)
{
	//Gauss 2
	if (N == 2)
	{
		double Gauss[2] = { -0.5773502691896250,0.5773502691896250 };
		double GaussWeight[2] = { 1.000000,1.0000000 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 3
	if (N == 3)
	{
		double Gauss[3] = { -0.774596669241484,0,0.774596669241484 };
		double GaussWeight[3] = { 5 / 9 ,8 / 9, 5 / 9 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 5
	if (N == 5)
	{
		double Gauss[5] = { -0.9061798459386640 ,-0.5384693101056830, 0.0000000000000000, 0.5384693101056830,0.9061798459386640 };
		double GaussWeight[5] = { 0.2369268850561890, 0.4786286704993660, 0.5688888888888880,0.4786286704993660,0.2369268850561890 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 8
	if (N == 8)
	{
		double Gauss[8] = { -0.960289856497536,	-0.796666477413627,	-0.525532409916329,	-0.183434642495650,	0.183434642495650,	0.525532409916329,	0.796666477413627,	0.960289856497536 };
		double GaussWeight[8] = { 0.101228536290376,	0.222381034453374,	0.313706645877887,	0.362683783378362,	0.362683783378362,	0.313706645877887,	0.222381034453374,	0.101228536290376 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 10
	if (N == 10)
	{
		double Gauss[10] = { -0.9739065285171710,-0.8650633666889840,-0.6794095682990240,-0.4333953941292470,-0.1488743389816310 ,0.148874338981631,0.433395394129247,0.679409568299024,0.865063366688984, 0.973906528517171 };
		double GaussWeight[10] = { 0.0666713443086881,0.1494513491505800,0.2190863625159820,0.2692667193099960,0.295524224714752,0.295524224714752,0.269266719309996,0.219086362515982,0.14945134915058,0.0666713443086881 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 15
	if (N == 15)
	{
		double Gauss[15] = { -0.987992518020486,	-0.937273392400706,	-0.848206583410427,	-0.724417731360170,	-0.570972172608539,	-0.394151347077563,	-0.201194093997435,	0,	0.201194093997435,	0.394151347077563,	0.570972172608539,	0.724417731360170,	0.848206583410427,	0.937273392400706,	0.987992518020486 };
		double GaussWeight[15] = { 0.0307532419961171,	0.0703660474881080,	0.107159220467172,	0.139570677926154,	0.166269205816994,	0.186161000015562,	0.198431485327112,	0.202578241925561,	0.198431485327112,	0.186161000015562,	0.166269205816994,	0.139570677926154,	0.107159220467172,	0.0703660474881080,	0.0307532419961171 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 20
	if (N == 20)
	{
		double Gauss[20];
		double GaussWeight[20];
		Gauss[0] = -0.9931285991850940; GaussWeight[0] = 0.0176140071391521; Gauss[1] = -0.9639719272779130; GaussWeight[1] = 0.0406014298003869; Gauss[2] = -0.9122344282513260; GaussWeight[2] = 0.0626720483341090; Gauss[3] = -0.8391169718222180; GaussWeight[3] = 0.0832767415767047; Gauss[4] = -0.746331906460150; GaussWeight[4] = 0.1019301198172400; Gauss[5] = -0.636053680726515; GaussWeight[5] = 0.1181945319615180; Gauss[6] = -0.5108670019508270; GaussWeight[6] = 0.1316886384491760; Gauss[7] = -0.3737060887154190; GaussWeight[7] = 0.1420961093183820; Gauss[8] = -0.227785851141645; GaussWeight[8] = 0.1491729864726030; Gauss[9] = -0.0765265211334973; GaussWeight[9] = 0.1527533871307250; Gauss[10] = 0.0765265211334973; GaussWeight[10] = 0.1527533871307250; Gauss[11] = 0.2277858511416450; GaussWeight[11] = 0.1491729864726030; Gauss[12] = 0.373706088715419; GaussWeight[12] = 0.1420961093183820; Gauss[13] = 0.5108670019508270; GaussWeight[13] = 0.1316886384491760; Gauss[14] = 0.6360536807265150; GaussWeight[14] = 0.1181945319615180; Gauss[15] = 0.7463319064601500; GaussWeight[15] = 0.1019301198172400; Gauss[16] = 0.8391169718222180; GaussWeight[16] = 0.0832767415767047; Gauss[17] = 0.9122344282513260; GaussWeight[17] = 0.0626720483341090; Gauss[18] = 0.9639719272779130; GaussWeight[18] = 0.0406014298003869; Gauss[19] = 0.9931285991850940; GaussWeight[19] = 0.0176140071391521;
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 25
	if (N == 25)
	{
		double Gauss[25] = { -0.995556969790498	,-0.976663921459518	,-0.942974571228974	,-0.894991997878275	,-0.833442628760834	,-0.759259263037358	,-0.673566368473468	,-0.577662930241223	,-0.473002731445715	,-0.361172305809388	,-0.243866883720988	,-0.122864692610710,	0,	0.122864692610710,	0.243866883720988,	0.361172305809388,	0.473002731445715	,0.577662930241223	,0.673566368473468	,0.759259263037358	,0.833442628760834	,0.894991997878275,	0.942974571228974,	0.976663921459518,	0.995556969790498 };
		double GaussWeight[25] = { 0.0113937985010262	,0.0263549866150321	,0.0409391567013064	,0.0549046959758352	,0.0680383338123569	,0.0801407003350011	,0.0910282619829636	,0.100535949067051	,0.108519624474264	,0.114858259145712	,0.119455763535785	,0.122242442990310	,0.123176053726715	,0.122242442990310	,0.119455763535785	,0.114858259145712	,0.108519624474264	,0.100535949067051	,0.0910282619829636	,0.0801407003350011	,0.0680383338123569	,0.0549046959758352	,0.0409391567013064	,0.0263549866150321	,0.0113937985010262 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 30
	if (N == 30)
	{
		double Gauss[30] = { -0.996893484074650	,-0.983668123279747	,-0.960021864968308	,-0.926200047429274	,-0.882560535792053	,-0.829565762382768	,-0.767777432104826	,-0.697850494793316	,-0.620526182989243	,-0.536624148142020	,-0.447033769538089	,-0.352704725530878	,-0.254636926167890	,-0.153869913608584	,-0.0514718425553177	,0.0514718425553177	,0.153869913608584	,0.254636926167890	,0.352704725530878	,0.447033769538089	,0.536624148142020	,0.620526182989243	,0.697850494793316	,0.767777432104826	,0.829565762382768	,0.882560535792053	,0.926200047429274	,0.960021864968308	,0.983668123279747	,0.996893484074650 };
		double GaussWeight[30] = { 0.00796819249616665	,0.0184664683110910	,0.0287847078833234	,0.0387991925696271	,0.0484026728305940	,0.0574931562176190	,0.0659742298821804	,0.0737559747377051	,0.0807558952294202	,0.0868997872010830	,0.0921225222377862	,0.0963687371746442	,0.0995934205867953	,0.101762389748405	,0.102852652893559	,0.102852652893559	,0.101762389748405	,0.0995934205867953	,0.0963687371746442	,0.0921225222377862	,0.0868997872010830	,0.0807558952294202	,0.0737559747377051	,0.0659742298821804	,0.0574931562176190	,0.0484026728305940	,0.0387991925696271	,0.0287847078833234	,0.0184664683110910	,0.00796819249616665 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 40
	if (N == 40)
	{
		double Gauss[40] = { -9.982377097105591e-01,-9.907262386994570e-01,-9.772599499837743e-01,-9.579168192137917e-01,-9.328128082786765e-01,-9.020988069688743e-01,-8.659595032122595e-01,-8.246122308333117e-01,-7.783056514265194e-01,-7.273182551899271e-01,-6.719566846141796e-01,-6.125538896679802e-01,-5.494671250951282e-01,-4.830758016861787e-01,-4.137792043716050e-01,-3.419940908257585e-01,-2.681521850072537e-01,-1.926975807013711e-01,-1.160840706752552e-01,-3.877241750605082e-02,3.877241750605082e-02,1.160840706752552e-01,1.926975807013711e-01,2.681521850072537e-01,3.419940908257585e-01,4.137792043716050e-01,4.830758016861787e-01,5.494671250951282e-01,6.125538896679802e-01,6.719566846141796e-01,7.273182551899271e-01,7.783056514265194e-01,8.246122308333117e-01,8.659595032122595e-01, 9.020988069688743e-01   ,  9.328128082786765e-01 ,9.579168192137917e-01   ,  9.772599499837743e-01   ,  9.907262386994570e-01   ,  9.982377097105591e-01 };
		double GaussWeight[40] = { 4.521277098533279e-03,1.049828453115287e-02   ,  1.642105838190785e-02   ,  2.224584919416696e-02   ,  2.793700698002344e-02  ,   3.346019528254776e-02 ,3.878216797447204e-02  ,   4.387090818567326e-02    , 4.869580763507216e-02  ,   5.322784698393681e-02  ,   5.743976909939157e-02  ,   6.130624249292890e-02 ,6.480401345660104e-02 ,    6.791204581523386e-02  ,   7.061164739128678e-02   ,  7.288658239580413e-02   ,  7.472316905796823e-02   ,  7.611036190062634e-02 , 7.703981816424797e-02  ,   7.750594797842486e-02   ,  7.750594797842486e-02   ,  7.703981816424797e-02  ,   7.611036190062634e-02   ,  7.472316905796823e-02 ,7.288658239580413e-02,     7.061164739128678e-02,     6.791204581523386e-02 ,    6.480401345660104e-02   ,  6.130624249292890e-02  ,   5.743976909939157e-02 ,5.322784698393681e-02 ,    4.869580763507216e-02  ,   4.387090818567326e-02 ,    3.878216797447204e-02 ,    3.346019528254776e-02   ,  2.793700698002344e-02 ,2.224584919416696e-02   ,  1.642105838190785e-02   ,  1.049828453115287e-02   ,  4.521277098533279e-03 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 50
	if (N == 50)
	{
		double Gauss[50] = { -0.998866404420071	,-0.994031969432091	,-0.985354084048006	,-0.972864385106692	,-0.956610955242808	,-0.936656618944878	,-0.913078556655792	,-0.885967979523613	,-0.855429769429946	,-0.821582070859336	,-0.784555832900399	,-0.744494302226068	,-0.701552468706822	,-0.655896465685439	,-0.607702927184950	,-0.557158304514650	,-0.504458144907464	,-0.449806334974039	,-0.393414311897565	,-0.335500245419437	,-0.276288193779532	,-0.216007236876042	,-0.154890589998146	,-0.0931747015600861	,-0.0310983383271889	,0.0310983383271889	,0.0931747015600861	,0.154890589998146	,0.216007236876042	,0.276288193779532	,0.335500245419437	,0.393414311897565	,0.449806334974039	,0.504458144907464	,0.557158304514650	,0.607702927184950	,0.655896465685439	,0.701552468706822	,0.744494302226068	,0.784555832900399	,0.821582070859336	,0.855429769429946	,0.885967979523613	,0.913078556655792	,0.936656618944878	,0.956610955242808	,0.972864385106692	,0.985354084048006	,0.994031969432091	,0.998866404420071 };
		double GaussWeight[50] = { 0.00290862255315528	,0.00675979919574537	,0.0105905483836509	,0.0143808227614856	,0.0181155607134894	,0.0217802431701248	,0.0253606735700124	,0.0288429935805353	,0.0322137282235780	,0.0354598356151462	,0.0385687566125877	,0.0415284630901477	,0.0443275043388032	,0.0469550513039484	,0.0494009384494664	,0.0516557030695811	,0.0537106218889963	,0.0555577448062126	,0.0571899256477284	,0.0586008498132225	,0.0597850587042655	,0.0607379708417702	,0.0614558995903167	,0.0619360674206832	,0.0621766166553473	,0.0621766166553473	,0.0619360674206832	,0.0614558995903167	,0.0607379708417702	,0.0597850587042655	,0.0586008498132225	,0.0571899256477284	,0.0555577448062126	,0.0537106218889963	,0.0516557030695811	,0.0494009384494664	,0.0469550513039484	,0.0443275043388032	,0.0415284630901477	,0.0385687566125877	,0.0354598356151462	,0.0322137282235780	,0.0288429935805353	,0.0253606735700124	,0.0217802431701248	,0.0181155607134894	,0.0143808227614856	,0.0105905483836509	,0.00675979919574537	,0.00290862255315528 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 60
	if (N == 60)
	{
		double Gauss[60] = { -0.999210123227436	,-0.995840525118838	,-0.989787895222222	,-0.981067201752598	,-0.969701788765053	,-0.955722255839996	,-0.939166276116423	,-0.920078476177628	,-0.898510310810046	,-0.874519922646898	,-0.848171984785930	,-0.819537526162146	,-0.788693739932264	,-0.755723775306586	,-0.720716513355730	,-0.683766327381356	,-0.644972828489477	,-0.604440597048510	,-0.562278900753945	,-0.518601400058570	,-0.473525841761707	,-0.427173741583078	,-0.379670056576798	,-0.331142848268448	,-0.281722937423262	,-0.231543551376029	,-0.180739964873425	,-0.129449135396945	,-0.0778093339495366	,-0.0259597723012478	,0.0259597723012478	,0.0778093339495366	,0.129449135396945	,0.180739964873425	,0.231543551376029	,0.281722937423262	,0.331142848268448	,0.379670056576798	,0.427173741583078	,0.473525841761707	,0.518601400058570	,0.562278900753945	,0.604440597048510	,0.644972828489477	,0.683766327381356	,0.720716513355730	,0.755723775306586	,0.788693739932264	,0.819537526162146	,0.848171984785930	,0.874519922646898	,0.898510310810046	,0.920078476177628	,0.939166276116423	,0.955722255839996	,0.969701788765053	,0.981067201752598	,0.989787895222222	,0.995840525118838	,0.999210123227436 };
		double GaussWeight[60] = { 0.00202681196887330	,0.00471272992695357	,0.00738993116334556	,0.0100475571822880	,0.0126781664768160	,0.0152746185967848	,0.0178299010142077	,0.0203371207294573	,0.0227895169439979	,0.0251804776215213	,0.0275035567499248	,0.0297524915007889	,0.0319212190192963	,0.0340038927249464	,0.0359948980510845	,0.0378888675692434	,0.0396806954523808	,0.0413655512355848	,0.0429388928359356	,0.0443964787957871	,0.0457343797161146	,0.0469489888489122	,0.0480370318199712	,0.0489955754557568	,0.0498220356905502	,0.0505141845325093	,0.0510701560698557	,0.0514884515009809	,0.0517679431749101	,0.0519078776312206	,0.0519078776312206	,0.0517679431749101	,0.0514884515009809	,0.0510701560698557	,0.0505141845325093	,0.0498220356905502	,0.0489955754557568	,0.0480370318199712	,0.0469489888489122	,0.0457343797161146	,0.0443964787957871	,0.0429388928359356	,0.0413655512355848	,0.0396806954523808	,0.0378888675692434	,0.0359948980510845	,0.0340038927249464	,0.0319212190192963	,0.0297524915007889	,0.0275035567499248	,0.0251804776215213	,0.0227895169439979	,0.0203371207294573	,0.0178299010142077	,0.0152746185967848	,0.0126781664768160	,0.0100475571822880	,0.00738993116334556	,0.00471272992695357	,0.00202681196887330 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 70
	if (N == 70)
	{
		double Gauss[70] = { -0.999418285973576	,-0.996936251961680	,-0.992476055211690	,-0.986045558070399	,-0.977657405957592	,-0.967328223664986	,-0.955078509114293	,-0.940932579003815	,-0.924918516897935	,-0.907068116260923	,-0.887416816863348	,-0.866003634213859	,-0.842871081998980	,-0.818065087625441	,-0.791634901007893	,-0.763632996771900	,-0.734114970060943	,-0.703139426151529	,-0.670767864094077	,-0.637064554609778	,-0.602096412485356	,-0.565932863718808	,-0.528645707679711	,-0.490308974557637	,-0.450998778381648	,-0.410793165902631	,-0.369771961638462	,-0.328016609389643	,-0.285610010540038	,-0.242636359463741	,-0.199180976364858	,-0.155330137882070	,-0.111170905794299	,-0.0667909541675513	,-0.0222783952861403	,0.0222783952861403	,0.0667909541675513	,0.111170905794299	,0.155330137882070	,0.199180976364858	,0.242636359463741	,0.285610010540038	,0.328016609389643	,0.369771961638462	,0.410793165902631	,0.450998778381648	,0.490308974557637	,0.528645707679711	,0.565932863718808	,0.602096412485356	,0.637064554609778	,0.670767864094077	,0.703139426151529	,0.734114970060943	,0.763632996771900	,0.791634901007893	,0.818065087625441	,0.842871081998980	,0.866003634213859	,0.887416816863348	,0.907068116260923	,0.924918516897935	,0.940932579003815	,0.955078509114293	,0.967328223664986	,0.977657405957592	,0.986045558070399	,0.992476055211690	,0.996936251961680	,0.999418285973576 };
		double GaussWeight[70] = { 0.00149272128884452	,0.00347189489307816	,0.00544711187421725	,0.00741176936319022	,0.00936176276969899	,0.0112931846499316	,0.0132021908146768	,0.0150849878654431	,0.0169378363763029	,0.0187570570931334	,0.0205390378243264	,0.0222802404522566	,0.0239772078891003	,0.0256265709084685	,0.0272250548186644	,0.0287694859558083	,0.0302567979801542	,0.0316840379613085	,0.0330483722393725	,0.0343470920499065	,0.0355776189012924	,0.0367375096936727	,0.0378244615692228	,0.0388363164840734	,0.0397710654927766	,0.0406268527367896	,0.0414019791290452	,0.0420949057272844	,0.0427042567894498	,0.0432288225050687	,0.0436675613972014	,0.0440196023901835	,0.0442842465390554	,0.0444609684172463	,0.0445494171597546	,0.0445494171597546	,0.0444609684172463	,0.0442842465390554	,0.0440196023901835	,0.0436675613972014	,0.0432288225050687	,0.0427042567894498	,0.0420949057272844	,0.0414019791290452	,0.0406268527367896	,0.0397710654927766	,0.0388363164840734	,0.0378244615692228	,0.0367375096936727	,0.0355776189012924	,0.0343470920499065	,0.0330483722393725	,0.0316840379613085	,0.0302567979801542	,0.0287694859558083	,0.0272250548186644	,0.0256265709084685	,0.0239772078891003	,0.0222802404522566	,0.0205390378243264	,0.0187570570931334	,0.0169378363763029	,0.0150849878654431	,0.0132021908146768	,0.0112931846499316	,0.00936176276969899	,0.00741176936319022	,0.00544711187421725	,0.00347189489307816	,0.00149272128884452 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 80
	if (N == 80)
	{
		double Gauss[80] = { -0.999553822651631	,-0.997649864398238	,-0.994227540965688	,-0.989291302499756	,-0.982848572738629	,-0.974909140585728	,-0.965485089043799	,-0.954590766343635	,-0.942242761309873	,-0.928459877172446	,-0.913263102571758	,-0.896675579438771	,-0.878722567678214	,-0.859431406663111	,-0.838831473580255	,-0.816954138681464	,-0.793832717504606	,-0.769502420135041	,-0.744000297583597	,-0.717365185362100	,-0.689637644342028	,-0.660859898986120	,-0.631075773046872	,-0.600330622829752	,-0.568671268122710	,-0.536145920897132	,-0.502804111888785	,-0.468696615170545	,-0.433875370831756	,-0.398393405881969	,-0.362304753499487	,-0.325664370747702	,-0.288528054884512	,-0.250952358392272	,-0.212994502857666	,-0.174712291832647	,-0.136164022809144	,-0.0974083984415846	,-0.0585044371524207	,-0.0195113832567940	,0.0195113832567940	,0.0585044371524207	,0.0974083984415846	,0.136164022809144	,0.174712291832647	,0.212994502857666	,0.250952358392272	,0.288528054884512	,0.325664370747702	,0.362304753499487	,0.398393405881969	,0.433875370831756	,0.468696615170545	,0.502804111888785	,0.536145920897132	,0.568671268122710	,0.600330622829752	,0.631075773046872	,0.660859898986120	,0.689637644342028	,0.717365185362100	,0.744000297583597	,0.769502420135041	,0.793832717504606	,0.816954138681464	,0.838831473580255	,0.859431406663111	,0.878722567678214	,0.896675579438771	,0.913263102571758	,0.928459877172446	,0.942242761309873	,0.954590766343635	,0.965485089043799	,0.974909140585728	,0.982848572738629	,0.989291302499756	,0.994227540965688	,0.997649864398238	,0.999553822651631 };
		double GaussWeight[80] = { 0.00114495000318693	,0.00266353358951265	,0.00418031312469485	,0.00569092245140326	,0.00719290476811730	,0.00868394526926085	,0.0101617660411030	,0.0116241141207978	,0.0130687615924013	,0.0144935080405090	,0.0158961835837257	,0.0172746520562693	,0.0186268142082990	,0.0199506108781419	,0.0212440261157820	,0.0225050902463325	,0.0237318828659301	,0.0249225357641155	,0.0260752357675651	,0.0271882275004864	,0.0282598160572768	,0.0292883695832679	,0.0302723217595580	,0.0312101741881147	,0.0321004986734878	,0.0329419393976454	,0.0337332149846114	,0.0344731204517540	,0.0351605290447476	,0.0357943939534160	,0.0363737499058360	,0.0368977146382761	,0.0373654902387304	,0.0377763643620014	,0.0381297113144776	,0.0384249930069594	,0.0386617597740764	,0.0388396510590520	,0.0389583959627695	,0.0390178136563066	,0.0390178136563066	,0.0389583959627695	,0.0388396510590520	,0.0386617597740764	,0.0384249930069594	,0.0381297113144776	,0.0377763643620014	,0.0373654902387304	,0.0368977146382761	,0.0363737499058360	,0.0357943939534160	,0.0351605290447476	,0.0344731204517540	,0.0337332149846114	,0.0329419393976454	,0.0321004986734878	,0.0312101741881147	,0.0302723217595580	,0.0292883695832679	,0.0282598160572768	,0.0271882275004864	,0.0260752357675651	,0.0249225357641155	,0.0237318828659301	,0.0225050902463325	,0.0212440261157820	,0.0199506108781419	,0.0186268142082990	,0.0172746520562693	,0.0158961835837257	,0.0144935080405090	,0.0130687615924013	,0.0116241141207978	,0.0101617660411030	,0.00868394526926085	,0.00719290476811730	,0.00569092245140326	,0.00418031312469485	,0.00266353358951265	,0.00114495000318693 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 90
	if (N == 90)
	{
		double Gauss[90] = { -0.999646971286639	,-0.998140379938568	,-0.995431812058345	,-0.991523928811063	,-0.986421365057833	,-0.980130251345148	,-0.972658162090193	,-0.964014098171506	,-0.954208473881500	,-0.943253103645358	,-0.931161187500432	,-0.917947295066586	,-0.903627347931303	,-0.888218600434746	,-0.871739618862903	,-0.854210259067072	,-0.835651642533377	,-0.816086130929481	,-0.795537299158248	,-0.774029906950334	,-0.751589869029638	,-0.728244223887390	,-0.704021101202391	,-0.678949687946597	,-0.653060193216842	,-0.626383811835045	,-0.598952686760742	,-0.570799870361221	,-0.541959284585914	,-0.512465680093028	,-0.482354594377666	,-0.451662308951869	,-0.420425805628198	,-0.388682721959498	,-0.356471305888568	,-0.323830369662346	,-0.290799243066167	,-0.257417726034420	,-0.223726040694723	,-0.189764782903379	,-0.155574873330529	,-0.121197508153924	,-0.0866741094207348	,-0.0520462751372070	,-0.0173557291462996	,0.0173557291462996	,0.0520462751372070	,0.0866741094207348	,0.121197508153924	,0.155574873330529	,0.189764782903379	,0.223726040694723	,0.257417726034420	,0.290799243066167	,0.323830369662346	,0.356471305888568	,0.388682721959498	,0.420425805628198	,0.451662308951869	,0.482354594377666	,0.512465680093028	,0.541959284585914	,0.570799870361221	,0.598952686760742	,0.626383811835045	,0.653060193216842	,0.678949687946597	,0.704021101202391	,0.728244223887390	,0.751589869029638	,0.774029906950334	,0.795537299158248	,0.816086130929481	,0.835651642533377	,0.854210259067072	,0.871739618862903	,0.888218600434746	,0.903627347931303	,0.917947295066586	,0.931161187500432	,0.943253103645358	,0.954208473881500	,0.964014098171506	,0.972658162090193	,0.980130251345148	,0.986421365057833	,0.991523928811063	,0.995431812058345	,0.998140379938568	,0.999646971286639 };
		double GaussWeight[90] = { 0.000905932371214690	,0.00210777877452633	,0.00330886724333602	,0.00450612361367503	,0.00569798156074735	,0.00688298320846328	,0.00805969494461999	,0.00922669695774197	,0.0103825823098933	,0.0115259578891480	,0.0126554458371681	,0.0137696851123372	,0.0148673330880433	,0.0159470671510067	,0.0170075862852227	,0.0180476126344603	,0.0190658930391373	,0.0200612005446396	,0.0210323358787225	,0.0219781288959341	,0.0228974399871632	,0.0237891614525287	,0.0246522188359049	,0.0254855722194432	,0.0262882174765146	,0.0270591874815480	,0.0277975532753023	,0.0285024251841615	,0.0291729538921007	,0.0298083314640313	,0.0304077923192869	,0.0309706141540809	,0.0314961188118187	,0.0319836731002186	,0.0324326895542556	,0.0328426271440075	,0.0332129919265514	,0.0335433376411243	,0.0338332662468317	,0.0340824284022539	,0.0342905238863750	,0.0344573019603243	,0.0345825616694969	,0.0346661520856882	,0.0347079724889500	,0.0347079724889500	,0.0346661520856882	,0.0345825616694969	,0.0344573019603243	,0.0342905238863750	,0.0340824284022539	,0.0338332662468317	,0.0335433376411243	,0.0332129919265514	,0.0328426271440075	,0.0324326895542556	,0.0319836731002186	,0.0314961188118187	,0.0309706141540809	,0.0304077923192869	,0.0298083314640313	,0.0291729538921007	,0.0285024251841615	,0.0277975532753023	,0.0270591874815480	,0.0262882174765146	,0.0254855722194432	,0.0246522188359049	,0.0237891614525287	,0.0228974399871632	,0.0219781288959341	,0.0210323358787225	,0.0200612005446396	,0.0190658930391373	,0.0180476126344603	,0.0170075862852227	,0.0159470671510067	,0.0148673330880433	,0.0137696851123372	,0.0126554458371681	,0.0115259578891480	,0.0103825823098933	,0.00922669695774197	,0.00805969494461999	,0.00688298320846328	,0.00569798156074735	,0.00450612361367503	,0.00330886724333602	,0.00210777877452633	,0.000905932371214690 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 100
	if (N == 100)
	{
		double Gauss[100] = { -0.999713726773441	,-0.998491950639596	,-0.996295134733125	,-0.993124937037443	,-0.988984395242992	,-0.983877540706057	,-0.977809358486918	,-0.970785775763707	,-0.962813654255816	,-0.953900782925492	,-0.944055870136256	,-0.933288535043080	,-0.921609298145334	,-0.909029570982530	,-0.895561644970727	,-0.881218679385018	,-0.866014688497165	,-0.849964527879591	,-0.833083879888401	,-0.815389238339176	,-0.796897892390315	,-0.777627909649496	,-0.757598118519707	,-0.736828089802021	,-0.715338117573056	,-0.693149199355802	,-0.670283015603141	,-0.646761908514129	,-0.622608860203708	,-0.597847470247179	,-0.572501932621381	,-0.546597012065094	,-0.520158019881763	,-0.493210789208191	,-0.465781649773358	,-0.437897402172032	,-0.409585291678302	,-0.380872981624630	,-0.351788526372422	,-0.322360343900529	,-0.292617188038472	,-0.262588120371504	,-0.232302481844974	,-0.201789864095736	,-0.171080080538603	,-0.140203137236114	,-0.109189203580061	,-0.0780685828134366	,-0.0468716824215916	,-0.0156289844215431	,0.0156289844215431	,0.0468716824215916	,0.0780685828134366	,0.109189203580061	,0.140203137236114	,0.171080080538603	,0.201789864095736	,0.232302481844974	,0.262588120371504	,0.292617188038472	,0.322360343900529	,0.351788526372422	,0.380872981624630	,0.409585291678302	,0.437897402172032	,0.465781649773358	,0.493210789208191	,0.520158019881763	,0.546597012065094	,0.572501932621381	,0.597847470247179	,0.622608860203708	,0.646761908514129	,0.670283015603141	,0.693149199355802	,0.715338117573056	,0.736828089802021	,0.757598118519707	,0.777627909649496	,0.796897892390315	,0.815389238339176	,0.833083879888401	,0.849964527879591	,0.866014688497165	,0.881218679385018	,0.895561644970727	,0.909029570982530	,0.921609298145334	,0.933288535043080	,0.944055870136256	,0.953900782925492	,0.962813654255816	,0.970785775763707	,0.977809358486918	,0.983877540706057	,0.988984395242992	,0.993124937037443	,0.996295134733125	,0.998491950639596	,0.999713726773441 };
		double GaussWeight[100] = { 0.000734634490505571	,0.00170939265351808	,0.00268392537155346	,0.00365596120132641	,0.00462445006342215	,0.00558842800386548	,0.00654694845084530	,0.00749907325546469	,0.00844387146966896	,0.00938041965369448	,0.0103078025748689	,0.0112251140231860	,0.0121314576629795	,0.0130259478929715	,0.0139077107037188	,0.0147758845274413	,0.0156296210775460	,0.0164680861761452	,0.0172904605683235	,0.0180959407221281	,0.0188837396133749	,0.0196530874944353	,0.0204032326462094	,0.0211334421125276	,0.0218430024162474	,0.0225312202563363	,0.0231974231852541	,0.0238409602659682	,0.0244612027079571	,0.0250575444815796	,0.0256294029102081	,0.0261762192395456	,0.0266974591835709	,0.0271926134465769	,0.0276611982207924	,0.0281027556591011	,0.0285168543223951	,0.0289030896011252	,0.0292610841106382	,0.0295904880599126	,0.0298909795933329	,0.0301622651051692	,0.0304040795264548	,0.0306161865839805	,0.0307983790311526	,0.0309504788504910	,0.0310723374275665	,0.0311638356962099	,0.0312248842548494	,0.0312554234538634	,0.0312554234538634	,0.0312248842548494	,0.0311638356962099	,0.0310723374275665	,0.0309504788504910	,0.0307983790311526	,0.0306161865839805	,0.0304040795264548	,0.0301622651051692	,0.0298909795933329	,0.0295904880599126	,0.0292610841106382	,0.0289030896011252	,0.0285168543223951	,0.0281027556591011	,0.0276611982207924	,0.0271926134465769	,0.0266974591835709	,0.0261762192395456	,0.0256294029102081	,0.0250575444815796	,0.0244612027079571	,0.0238409602659682	,0.0231974231852541	,0.0225312202563363	,0.0218430024162474	,0.0211334421125276	,0.0204032326462094	,0.0196530874944353	,0.0188837396133749	,0.0180959407221281	,0.0172904605683235	,0.0164680861761452	,0.0156296210775460	,0.0147758845274413	,0.0139077107037188	,0.0130259478929715	,0.0121314576629795	,0.0112251140231860	,0.0103078025748689	,0.00938041965369448	,0.00844387146966896	,0.00749907325546469	,0.00654694845084530	,0.00558842800386548	,0.00462445006342215	,0.00365596120132641	,0.00268392537155346	,0.00170939265351808	,0.000734634490505571 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 110
	if (N == 110)
	{
		double Gauss[110] = { -0.999763194109853	,-0.998752494089402	,-0.996935024690191	,-0.994311850823600	,-0.990885043159995	,-0.986657359588479	,-0.981632213192041	,-0.975813663934759	,-0.969206413903879	,-0.961815803018185	,-0.953647804520909	,-0.944709020068334	,-0.935006674354058	,-0.924548609249093	,-0.913343277452476	,-0.901399735652954	,-0.888727637205083	,-0.875337224324529	,-0.861239319808274	,-0.846445318286095	,-0.830967177010178	,-0.814817406190207	,-0.798009058881658	,-0.780555720435431	,-0.762471497517321	,-0.743771006706173	,-0.724469362679938	,-0.704582165999162	,-0.684125490497768	,-0.663115870291337	,-0.641570286413370	,-0.619506153090336	,-0.596941303666598	,-0.573893976190585	,-0.550382798673877	,-0.526426774035089	,-0.502045264740751	,-0.477257977155578	,-0.452084945614786	,-0.426546516231330	,-0.400663330451145	,-0.374456308369683	,-0.347946631823240	,-0.321155727268719	,-0.294105248465692	,-0.266817058974733	,-0.239313214486187	,-0.211615944993647	,-0.183747636826553	,-0.155730814556430	,-0.127588122791390	,-0.0993423078736232	,-0.0710161994946515	,-0.0426326922432233	,-0.0142147271007502	,0.0142147271007502	,0.0426326922432233	,0.0710161994946515	,0.0993423078736232	,0.127588122791390	,0.155730814556430	,0.183747636826553	,0.211615944993647	,0.239313214486187	,0.266817058974733	,0.294105248465692	,0.321155727268719	,0.347946631823240	,0.374456308369683	,0.400663330451145	,0.426546516231330	,0.452084945614786	,0.477257977155578	,0.502045264740751	,0.526426774035089	,0.550382798673877	,0.573893976190585	,0.596941303666598	,0.619506153090336	,0.641570286413370	,0.663115870291337	,0.684125490497768	,0.704582165999162	,0.724469362679938	,0.743771006706173	,0.762471497517321	,0.780555720435431	,0.798009058881658	,0.814817406190207	,0.830967177010178	,0.846445318286095	,0.861239319808274	,0.875337224324529	,0.888727637205083	,0.901399735652954	,0.913343277452476	,0.924548609249093	,0.935006674354058	,0.944709020068334	,0.953647804520909	,0.961815803018185	,0.969206413903879	,0.975813663934759	,0.981632213192041	,0.986657359588479	,0.990885043159995	,0.994311850823600	,0.996935024690191	,0.998752494089402	,0.999763194109853 };
		double GaussWeight[110] = { 0.000607696399185895	,0.00141412487600019	,0.00222060546191683	,0.00302539507433425	,0.00382776026499883	,0.00462703791705988	,0.00542257814058416	,0.00621373662242899	,0.00699987337810134	,0.00778035277100900	,0.00855454385793997	,0.00932182083395445	,0.0100815635102698	,0.0108331578025071	,0.0115759962205267	,0.0123094783559997	,0.0130330113657858	,0.0137460104500006	,0.0144478993240229	,0.0151381106838707	,0.0158160866644674	,0.0164812792903681	,0.0171331509185498	,0.0177711746728855	,0.0183948348699397	,0.0190036274357330	,0.0195970603131338	,0.0201746538595445	,0.0207359412345593	,0.0212804687772773	,0.0218077963729667	,0.0223174978087821	,0.0228091611182459	,0.0232823889142181	,0.0237367987100806	,0.0241720232288814	,0.0245877107001836	,0.0249835251443833	,0.0253591466442638	,0.0257142716035689	,0.0260486129923825	,0.0263619005791213	,0.0266538811489487	,0.0269243187084356	,0.0271729946763015	,0.0273997080600830	,0.0276042756185859	,0.0277865320099894	,0.0279463299254841	,0.0280835402083338	,0.0281980519582665	,0.0282897726211096	,0.0283586280635972	,0.0284045626332874	,0.0284275392035439	,0.0284275392035439	,0.0284045626332874	,0.0283586280635972	,0.0282897726211096	,0.0281980519582665	,0.0280835402083338	,0.0279463299254841	,0.0277865320099894	,0.0276042756185859	,0.0273997080600830	,0.0271729946763015	,0.0269243187084356	,0.0266538811489487	,0.0263619005791213	,0.0260486129923825	,0.0257142716035689	,0.0253591466442638	,0.0249835251443833	,0.0245877107001836	,0.0241720232288814	,0.0237367987100806	,0.0232823889142181	,0.0228091611182459	,0.0223174978087821	,0.0218077963729667	,0.0212804687772773	,0.0207359412345593	,0.0201746538595445	,0.0195970603131338	,0.0190036274357330	,0.0183948348699397	,0.0177711746728855	,0.0171331509185498	,0.0164812792903681	,0.0158160866644674	,0.0151381106838707	,0.0144478993240229	,0.0137460104500006	,0.0130330113657858	,0.0123094783559997	,0.0115759962205267	,0.0108331578025071	,0.0100815635102698	,0.00932182083395445	,0.00855454385793997	,0.00778035277100900	,0.00699987337810134	,0.00621373662242899	,0.00542257814058416	,0.00462703791705988	,0.00382776026499883	,0.00302539507433425	,0.00222060546191683	,0.00141412487600019	,0.000607696399185895 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 120
	if (N == 120)
	{
		double Gauss[120] = { -0.999800865658959	,-0.998950921674035	,-0.997422413613291	,-0.995216040598908	,-0.992333260616166	,-0.988776023071526	,-0.984546742413411	,-0.979648291821042	,-0.974084000010228	,-0.967857648557972	,-0.960973469171583	,-0.953436140742276	,-0.945250786131476	,-0.936422968671676	,-0.926958688375774	,-0.916864377853624	,-0.906146897936738	,-0.894813533013157	,-0.882871986075170	,-0.870330373482944	,-0.857197219447483	,-0.843481450236549	,-0.829192388107437	,-0.814339744970676	,-0.798933615788958	,-0.782984471715736	,-0.766503152978177	,-0.749500861509276	,-0.731989153334141	,-0.713979930715622	,-0.695485434064622	,-0.676518233620576	,-0.657091220907766	,-0.637217599973260	,-0.616910878412454	,-0.596184858188287	,-0.575053626250394	,-0.553531544960556	,-0.531633242330966	,-0.509373602081932	,-0.486767753525790	,-0.463831061283892	,-0.440579114843658	,-0.417027717962798	,-0.393192877927895	,-0.369090794674659	,-0.344737849777241	,-0.320150595314089	,-0.295345742617922	,-0.270340150917462	,-0.245150815878639	,-0.219794858053070	,-0.194289511241657	,-0.168652110781204	,-0.142900081762031	,-0.117050927184584	,-0.0911222160630847	,-0.0651315714843215	,-0.0390966586296822	,-0.0130351727685786	,0.0130351727685786	,0.0390966586296822	,0.0651315714843215	,0.0911222160630847	,0.117050927184584	,0.142900081762031	,0.168652110781204	,0.194289511241657	,0.219794858053070	,0.245150815878639	,0.270340150917462	,0.295345742617922	,0.320150595314089	,0.344737849777241	,0.369090794674659	,0.393192877927895	,0.417027717962798	,0.440579114843658	,0.463831061283892	,0.486767753525790	,0.509373602081932	,0.531633242330966	,0.553531544960556	,0.575053626250394	,0.596184858188287	,0.616910878412454	,0.637217599973260	,0.657091220907766	,0.676518233620576	,0.695485434064622	,0.713979930715622	,0.731989153334141	,0.749500861509276	,0.766503152978177	,0.782984471715736	,0.798933615788958	,0.814339744970676	,0.829192388107437	,0.843481450236549	,0.857197219447483	,0.870330373482944	,0.882871986075170	,0.894813533013157	,0.906146897936738	,0.916864377853624	,0.926958688375774	,0.936422968671676	,0.945250786131476	,0.953436140742276	,0.960973469171583	,0.967857648557972	,0.974084000010228	,0.979648291821042	,0.984546742413411	,0.988776023071526	,0.992333260616166	,0.995216040598908	,0.997422413613291	,0.998950921674035	,0.999800865658959 };
		double GaussWeight[120] = { 0.000511026063694716	,0.00118923432779327	,0.00186763923077199	,0.00254486205401025	,0.00322037273392629	,0.00389369986290145	,0.00456438254088213	,0.00523196381429558	,0.00589598949928373	,0.00655600807167165	,0.00721157083237960	,0.00786223215769001	,0.00850754977865233	,0.00914708507073764	,0.00978040334649752	,0.0104070741481252	,0.0110266715384281	,0.0116387743894037	,0.0122429666679165	,0.0128388377181243	,0.0134259825403743	,0.0140040020663295	,0.0145725034301117	,0.0151311002352592	,0.0156794128173061	,0.0162170685018000	,0.0167437018575773	,0.0172589549451213	,0.0177624775598333	,0.0182539274700487	,0.0187329706496370	,0.0191992815050258	,0.0196525430964950	,0.0200924473535893	,0.0205186952845038	,0.0209309971792998	,0.0213290728068114	,0.0217126516051107	,0.0220814728654004	,0.0224352859092113	,0.0227738502587807	,0.0230969358004986	,0.0234043229413100	,0.0236958027579665	,0.0239711771390253	,0.0242302589195009	,0.0244728720080758	,0.0246988515067856	,0.0249080438230954	,0.0251003067742927	,0.0252755096841254	,0.0254335334716194	,0.0255742707320141	,0.0256976258097633	,0.0258035148635495	,0.0258918659232686	,0.0259626189389465	,0.0260157258215529	,0.0260511504756872	,0.0260688688241102	,0.0260688688241102	,0.0260511504756872	,0.0260157258215529	,0.0259626189389465	,0.0258918659232686	,0.0258035148635495	,0.0256976258097633	,0.0255742707320141	,0.0254335334716194	,0.0252755096841254	,0.0251003067742927	,0.0249080438230954	,0.0246988515067856	,0.0244728720080758	,0.0242302589195009	,0.0239711771390253	,0.0236958027579665	,0.0234043229413100	,0.0230969358004986	,0.0227738502587807	,0.0224352859092113	,0.0220814728654004	,0.0217126516051107	,0.0213290728068114	,0.0209309971792998	,0.0205186952845038	,0.0200924473535893	,0.0196525430964950	,0.0191992815050258	,0.0187329706496370	,0.0182539274700487	,0.0177624775598333	,0.0172589549451213	,0.0167437018575773	,0.0162170685018000	,0.0156794128173061	,0.0151311002352592	,0.0145725034301117	,0.0140040020663295	,0.0134259825403743	,0.0128388377181243	,0.0122429666679165	,0.0116387743894037	,0.0110266715384281	,0.0104070741481252	,0.00978040334649752	,0.00914708507073764	,0.00850754977865233	,0.00786223215769001	,0.00721157083237960	,0.00655600807167165	,0.00589598949928373	,0.00523196381429558	,0.00456438254088213	,0.00389369986290145	,0.00322037273392629	,0.00254486205401025	,0.00186763923077199	,0.00118923432779327	,0.000511026063694716 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 160
	if (N == 160)
	{
		double Gauss[160] = { -0.999887752272163	,-0.999408620641306	,-0.998546818751517	,-0.997302485277280	,-0.995676073488025	,-0.993668200892097	,-0.991279634879646	,-0.988511289796807	,-0.985364225895890	,-0.981839648696526	,-0.977938908433306	,-0.973663499498934	,-0.969015059852788	,-0.963995370383663	,-0.958606354222183	,-0.952850076000975	,-0.946728741061852	,-0.940244694609783	,-0.933400420813678	,-0.926198541854185	,-0.918641816918757	,-0.910733141144340	,-0.902475544508024	,-0.893872190666099	,-0.884926375741905	,-0.875641527062964	,-0.866021201847856	,-0.856069085843345	,-0.845788991912262	,-0.835184858572704	,-0.824260748489090	,-0.813020846915656	,-0.801469460092987	,-0.789611013598204	,-0.777450050649423	,-0.764991230365148	,-0.752239325979266	,-0.739199223012315	,-0.725875917399735	,-0.712274513577825	,-0.698400222528118	,-0.684258359780950	,-0.669854343378962	,-0.655193691801340	,-0.640282021849564	,-0.625125046495499	,-0.609728572692629	,-0.594098499151298	,-0.578240814078786	,-0.562161592885103	,-0.545866995855374	,-0.529363265789701	,-0.512656725611419	,-0.495753775944651	,-0.478660892662094	,-0.461384624403979	,-0.443931590069149	,-0.426308476279219	,-0.408522034816796	,-0.390579080038724	,-0.372486486265367	,-0.354251185146911	,-0.335880163007704	,-0.317380458169650	,-0.298759158255680	,-0.280023397474338	,-0.261180353886511	,-0.242237246655374	,-0.223201333280566	,-0.204079906817702	,-0.184880293084242	,-0.165609847852822	,-0.146275954033103	,-0.126886018843225	,-0.107447470971946	,-0.0879677577325559	,-0.0684543422096580	,-0.0489147003999024	,-0.0293563183477774	,-0.00978668927754903	,0.00978668927754903	,0.0293563183477774	,0.0489147003999024	,0.0684543422096580	,0.0879677577325559	,0.107447470971946	,0.126886018843225	,0.146275954033103	,0.165609847852822	,0.184880293084242	,0.204079906817702	,0.223201333280566	,0.242237246655374	,0.261180353886511	,0.280023397474338	,0.298759158255680	,0.317380458169650	,0.335880163007704	,0.354251185146911	,0.372486486265367	,0.390579080038724	,0.408522034816796	,0.426308476279219	,0.443931590069149	,0.461384624403979	,0.478660892662094	,0.495753775944651	,0.512656725611419	,0.529363265789701	,0.545866995855374	,0.562161592885103	,0.578240814078786	,0.594098499151298	,0.609728572692629	,0.625125046495499	,0.640282021849564	,0.655193691801340	,0.669854343378962	,0.684258359780950	,0.698400222528118	,0.712274513577825	,0.725875917399735	,0.739199223012315	,0.752239325979266	,0.764991230365148	,0.777450050649423	,0.789611013598204	,0.801469460092987	,0.813020846915656	,0.824260748489090	,0.835184858572704	,0.845788991912262	,0.856069085843345	,0.866021201847856	,0.875641527062964	,0.884926375741905	,0.893872190666099	,0.902475544508024	,0.910733141144340	,0.918641816918757	,0.926198541854185	,0.933400420813678	,0.940244694609783	,0.946728741061852	,0.952850076000975	,0.958606354222183	,0.963995370383663	,0.969015059852788	,0.973663499498934	,0.977938908433306	,0.981839648696526	,0.985364225895890	,0.988511289796807	,0.991279634879646	,0.993668200892097	,0.995676073488025	,0.997302485277280	,0.998546818751517	,0.999408620641306	,0.999887752272163 };
		double GaussWeight[160] = { 0.000288058528521182	,0.000670438320152384	,0.00105312767221372	,0.00143546275348318	,0.00181725775966521	,0.00219835949679650	,0.00257862012183753	,0.00295789332503937	,0.00333603354782084	,0.00371289580258752	,0.00408833564856797	,0.00446220921633500	,0.00483437324970280	,0.00520468515444728	,0.00557300304986092	,0.00593918582148590	,0.00630309317427814	,0.00666458568583780	,0.00702352485951407	,0.00737977317727395	,0.00773319415226699	,0.00808365238103638	,0.00843101359534256	,0.00877514471356825	,0.00911591389168001	,0.00945319057372323	,0.00978684554182840	,0.0101167509657079	,0.0104427804516248	,0.0107648090908107	,0.0110827135073188	,0.0113963719052876	,0.0117056641156023	,0.0120104716419322	,0.0123106777061278	,0.0126061672929599	,0.0128968271941836	,0.0131825460519103	,0.0134632144012705	,0.0137387247123518	,0.0140089714313955	,0.0142738510212363	,0.0145332620009690	,0.0147871049848281	,0.0150352827202646	,0.0152777001252044	,0.0155142643244775	,0.0157448846853992	,0.0159694728524940	,0.0161879427813458	,0.0164002107715639	,0.0166061954988498	,0.0168058180461545	,0.0169990019339127	,0.0171856731493448	,0.0173657601748115	,0.0175391940152144	,0.0177059082244292	,0.0178658389307628	,0.0180189248614234	,0.0181651073659960	,0.0183043304389119	,0.0184365407409065	,0.0185616876194535	,0.0186797231281719	,0.0187906020451952	,0.0188942818904959	,0.0189907229421616	,0.0190798882516128	,0.0191617436577585	,0.0192362578000844	,0.0193034021306675	,0.0193631509251134	,0.0194154812924122	,0.0194603731837082	,0.0194978093999806	,0.0195277755986340	,0.0195502602989922	,0.0195652548866974	,0.0195727536170100	,0.0195727536170100	,0.0195652548866974	,0.0195502602989922	,0.0195277755986340	,0.0194978093999806	,0.0194603731837082	,0.0194154812924122	,0.0193631509251134	,0.0193034021306675	,0.0192362578000844	,0.0191617436577585	,0.0190798882516128	,0.0189907229421616	,0.0188942818904959	,0.0187906020451952	,0.0186797231281719	,0.0185616876194535	,0.0184365407409065	,0.0183043304389119	,0.0181651073659960	,0.0180189248614234	,0.0178658389307628	,0.0177059082244292	,0.0175391940152144	,0.0173657601748115	,0.0171856731493448	,0.0169990019339127	,0.0168058180461545	,0.0166061954988498	,0.0164002107715639	,0.0161879427813458	,0.0159694728524940	,0.0157448846853992	,0.0155142643244775	,0.0152777001252044	,0.0150352827202646	,0.0147871049848281	,0.0145332620009690	,0.0142738510212363	,0.0140089714313955	,0.0137387247123518	,0.0134632144012705	,0.0131825460519103	,0.0128968271941836	,0.0126061672929599	,0.0123106777061278	,0.0120104716419322	,0.0117056641156023	,0.0113963719052876	,0.0110827135073188	,0.0107648090908107	,0.0104427804516248	,0.0101167509657079	,0.00978684554182840	,0.00945319057372323	,0.00911591389168001	,0.00877514471356825	,0.00843101359534256	,0.00808365238103638	,0.00773319415226699	,0.00737977317727395	,0.00702352485951407	,0.00666458568583780	,0.00630309317427814	,0.00593918582148590	,0.00557300304986092	,0.00520468515444728	,0.00483437324970280	,0.00446220921633500	,0.00408833564856797	,0.00371289580258752	,0.00333603354782084	,0.00295789332503937	,0.00257862012183753	,0.00219835949679650	,0.00181725775966521	,0.00143546275348318	,0.00105312767221372	,0.000670438320152384	,0.000288058528521182 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 240
	if (N == 240)
	{
		double Gauss[240] = { -0.999950007741626	,-0.999736603137650	,-0.999352710098564	,-0.998798308846185	,-0.998073483532712	,-0.997178355342267	,-0.996113076181942	,-0.994877827483928	,-0.993472819863101	,-0.991898292977500	,-0.990154515447143	,-0.988241784790596	,-0.986160427365786	,-0.983910798309993	,-0.981493281476892	,-0.978908289369709	,-0.976156263070027	,-0.973237672162026	,-0.970153014652036	,-0.966902816883355	,-0.963487633446292	,-0.959908047083447	,-0.956164668590200	,-0.952258136710449	,-0.948189118027575	,-0.943958306850687	,-0.939566425096119	,-0.935014222164241	,-0.930302474811573	,-0.925431987018234	,-0.920403589850755	,-0.915218141320261	,-0.909876526236065	,-0.904379656054686	,-0.898728468724322	,-0.892923928524803	,-0.886967025903050	,-0.880858777304072	,-0.874600224997527	,-0.868192436899870	,-0.861636506392136	,-0.854933552133370	,-0.848084717869746	,-0.841091172239404	,-0.833954108573046	,-0.826674744690309	,-0.819254322691966	,-0.811694108747983	,-0.803995392881466	,-0.796159488748545	,-0.788187733414213	,-0.780081487124185	,-0.771842133072793	,-0.763471077166964	,-0.754969747786330	,-0.746339595539496	,-0.737582093016520	,-0.728698734537639	,-0.719691035898288	,-0.710560534110455	,-0.701308787140419	,-0.691937373642906	,-0.682447892691716	,-0.672841963506874	,-0.663121225178334	,-0.653287336386294	,-0.643341975118174	,-0.633286838382293	,-0.623123641918305	,-0.612854119904436	,-0.602480024661574	,-0.592003126354268	,-0.581425212688677	,-0.570748088607531	,-0.559973575982144	,-0.549103513301543	,-0.538139755358766	,-0.527084172934363	,-0.515938652477191	,-0.504705095782513	,-0.493385419667492	,-0.481981555644123	,-0.470495449589644	,-0.458929061414514	,-0.447284364727982	,-0.435563346501324	,-0.423768006728801	,-0.411900358086395	,-0.399962425588377	,-0.387956246241775	,-0.375883868698789	,-0.363747352907228	,-0.351548769759014	,-0.339290200736813	,-0.326973737558875	,-0.314601481822108	,-0.302175544643482	,-0.289698046299800	,-0.277171115865909	,-0.264596890851410	,-0.251977516835927	,-0.239315147102999	,-0.226611942272662	,-0.213870069932772	,-0.201091704269143	,-0.188279025694562	,-0.175434220476735	,-0.162559480365237	,-0.149657002217531	,-0.136728987624106	,-0.123777642532817	,-0.110805176872474	,-0.0978138041757565	,-0.0848057412015125	,-0.0717832075565031	,-0.0587484253166660	,-0.0457036186479554	,-0.0326510134268252	,-0.0195928368604209	,-0.00653131710654431	,0.00653131710654431	,0.0195928368604209	,0.0326510134268252	,0.0457036186479554	,0.0587484253166660	,0.0717832075565031	,0.0848057412015125	,0.0978138041757565	,0.110805176872474	,0.123777642532817	,0.136728987624106	,0.149657002217531	,0.162559480365237	,0.175434220476735	,0.188279025694562	,0.201091704269143	,0.213870069932772	,0.226611942272662	,0.239315147102999	,0.251977516835927	,0.264596890851410	,0.277171115865909	,0.289698046299800	,0.302175544643482	,0.314601481822108	,0.326973737558875	,0.339290200736813	,0.351548769759014	,0.363747352907228	,0.375883868698789	,0.387956246241775	,0.399962425588377	,0.411900358086395	,0.423768006728801	,0.435563346501324	,0.447284364727982	,0.458929061414514	,0.470495449589644	,0.481981555644123	,0.493385419667492	,0.504705095782513	,0.515938652477191	,0.527084172934363	,0.538139755358766	,0.549103513301543	,0.559973575982144	,0.570748088607531	,0.581425212688677	,0.592003126354268	,0.602480024661574	,0.612854119904436	,0.623123641918305	,0.633286838382293	,0.643341975118174	,0.653287336386294	,0.663121225178334	,0.672841963506874	,0.682447892691716	,0.691937373642906	,0.701308787140419	,0.710560534110455	,0.719691035898288	,0.728698734537639	,0.737582093016520	,0.746339595539496	,0.754969747786330	,0.763471077166964	,0.771842133072793	,0.780081487124185	,0.788187733414213	,0.796159488748545	,0.803995392881466	,0.811694108747983	,0.819254322691966	,0.826674744690309	,0.833954108573046	,0.841091172239404	,0.848084717869746	,0.854933552133370	,0.861636506392136	,0.868192436899870	,0.874600224997527	,0.880858777304072	,0.886967025903050	,0.892923928524803	,0.898728468724322	,0.904379656054686	,0.909876526236065	,0.915218141320261	,0.920403589850755	,0.925431987018234	,0.930302474811573	,0.935014222164241	,0.939566425096119	,0.943958306850687	,0.948189118027575	,0.952258136710449	,0.956164668590200	,0.959908047083447	,0.963487633446292	,0.966902816883355	,0.970153014652036	,0.973237672162026	,0.976156263070027	,0.978908289369709	,0.981493281476892	,0.983910798309993	,0.986160427365786	,0.988241784790596	,0.990154515447143	,0.991898292977500	,0.993472819863101	,0.994877827483928	,0.996113076181942	,0.997178355342267	,0.998073483532712	,0.998798308846185	,0.999352710098564	,0.999736603137650	,0.999950007741626 };
		double GaussWeight[240] = { 0.000128295209986303	,0.000298625582909719	,0.000469157256698254	,0.000639630788269057	,0.000809999577139253	,0.000980231471159985	,0.00115029660608832	,0.00132016568538638	,0.00148980961198420	,0.00165919938818522	,0.00182830608512461	,0.00199710083398567	,0.00216555482504024	,0.00233363930982096	,0.00250132560465756	,0.00266858509484865	,0.00283538923914393	,0.00300170957438313	,0.00316751772021437	,0.00333278538385168	,0.00349748436484845	,0.00366158655987407	,0.00382506396748581	,0.00398788869289108	,0.00415003295269598	,0.00431146907963863	,0.00447216952730438	,0.00463210687482230	,0.00479125383154115	,0.00494958324168387	,0.00510706808897974	,0.00526368150127294	,0.00541939675510697	,0.00557418728028377	,0.00572802666439689	,0.00588088865733803	,0.00603274717577571	,0.00618357630760570	,0.00633335031637235	,0.00648204364565971	,0.00662963092345234	,0.00677608696646453	,0.00692138678443710	,0.00706550558440164	,0.00720841877491101	,0.00735010197023521	,0.00749053099452241	,0.00762968188592442	,0.00776753090068469	,0.00790405451719049	,0.00803922943998589	,0.00817303260374716	,0.00830544117721796	,0.00843643256710560	,0.00856598442193585	,0.00869407463586711	,0.00882068135246201	,0.00894578296841722	,0.00906935813724954	,0.00919138577293814	,0.00931184505352282	,0.00943071542465675	,0.00954797660311374	,0.00966360858024932	,0.00977759162541467	,0.00988990628932340	,0.0100005334073706	,0.0101094541029024	,0.0102166497904371	,0.0103221021788369	,0.0104257932744283	,0.0105277053840729	,0.0106278211181863	,0.0107261233937054	,0.0108225954370030	,0.0109172207867504	,0.0110099832967261	,0.0111008671385707	,0.0111898568044878	,0.0112769371098905	,0.0113620931959917	,0.0114453105323402	,0.0115265749192997	,0.0116058724904719	,0.0116831897150623	,0.0117585134001893	,0.0118318306931352	,0.0119031290835392	,0.0119723964055324	,0.0120396208398134	,0.0121047909156652	,0.0121678955129125	,0.0122289238638190	,0.0122878655549248	,0.0123447105288237	,0.0123994490858784	,0.0124520718858768	,0.0125025699496245	,0.0125509346604778	,0.0125971577658139	,0.0126412313784386	,0.0126831479779325	,0.0127229004119342	,0.0127604818973607	,0.0127958860215647	,0.0128291067434289	,0.0128601383943969	,0.0128889756794404	,0.0129156136779626	,0.0129400478446379	,0.0129622740101873	,0.0129822883820903	,0.0130000875452315	,0.0130156684624838	,0.0130290284752259	,0.0130401653037971	,0.0130490770478848	,0.0130557621868505	,0.0130602195799876	,0.0130624484667171	,0.0130624484667171	,0.0130602195799876	,0.0130557621868505	,0.0130490770478848	,0.0130401653037971	,0.0130290284752259	,0.0130156684624838	,0.0130000875452315	,0.0129822883820903	,0.0129622740101873	,0.0129400478446379	,0.0129156136779626	,0.0128889756794404	,0.0128601383943969	,0.0128291067434289	,0.0127958860215647	,0.0127604818973607	,0.0127229004119342	,0.0126831479779325	,0.0126412313784386	,0.0125971577658139	,0.0125509346604778	,0.0125025699496245	,0.0124520718858768	,0.0123994490858784	,0.0123447105288237	,0.0122878655549248	,0.0122289238638190	,0.0121678955129125	,0.0121047909156652	,0.0120396208398134	,0.0119723964055324	,0.0119031290835392	,0.0118318306931352	,0.0117585134001893	,0.0116831897150623	,0.0116058724904719	,0.0115265749192997	,0.0114453105323402	,0.0113620931959917	,0.0112769371098905	,0.0111898568044878	,0.0111008671385707	,0.0110099832967261	,0.0109172207867504	,0.0108225954370030	,0.0107261233937054	,0.0106278211181863	,0.0105277053840729	,0.0104257932744283	,0.0103221021788369	,0.0102166497904371	,0.0101094541029024	,0.0100005334073706	,0.00988990628932340	,0.00977759162541467	,0.00966360858024932	,0.00954797660311374	,0.00943071542465675	,0.00931184505352282	,0.00919138577293814	,0.00906935813724954	,0.00894578296841722	,0.00882068135246201	,0.00869407463586711	,0.00856598442193585	,0.00843643256710560	,0.00830544117721796	,0.00817303260374716	,0.00803922943998589	,0.00790405451719049	,0.00776753090068469	,0.00762968188592442	,0.00749053099452241	,0.00735010197023521	,0.00720841877491101	,0.00706550558440164	,0.00692138678443710	,0.00677608696646453	,0.00662963092345234	,0.00648204364565971	,0.00633335031637235	,0.00618357630760570	,0.00603274717577571	,0.00588088865733803	,0.00572802666439689	,0.00557418728028377	,0.00541939675510697	,0.00526368150127294	,0.00510706808897974	,0.00494958324168387	,0.00479125383154115	,0.00463210687482230	,0.00447216952730438	,0.00431146907963863	,0.00415003295269598	,0.00398788869289108	,0.00382506396748581	,0.00366158655987407	,0.00349748436484845	,0.00333278538385168	,0.00316751772021437	,0.00300170957438313	,0.00283538923914393	,0.00266858509484865	,0.00250132560465756	,0.00233363930982096	,0.00216555482504024	,0.00199710083398567	,0.00182830608512461	,0.00165919938818522	,0.00148980961198420	,0.00132016568538638	,0.00115029660608832	,0.000980231471159985	,0.000809999577139253	,0.000639630788269057	,0.000469157256698254	,0.000298625582909719	,0.000128295209986303 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 320
	if (N == 320)
	{
		double Gauss[320] = { -0.999971849980221	,-0.999851682193268	,-0.999635503240465	,-0.999323285871991	,-0.998915054203942	,-0.998410846055053	,-0.997810709402343	,-0.997114701716636	,-0.996322889781677	,-0.995435349629393	,-0.994452166509907	,-0.993373434873440	,-0.992199258356498	,-0.990929749769477	,-0.989565031084503	,-0.988105233422948	,-0.986550497042385	,-0.984900971322824	,-0.983156814752184	,-0.981318194910950	,-0.979385288455988	,-0.977358281103520	,-0.975237367611243	,-0.973022751759588	,-0.970714646332118	,-0.968313273095077	,-0.965818862776063	,-0.963231655041857	,-0.960551898475386	,-0.957779850551837	,-0.954915777613911	,-0.951959954846234	,-0.948912666248912	,-0.945774204610239	,-0.942544871478572	,-0.939224977133349	,-0.935814840555278	,-0.932314789395691	,-0.928725159945059	,-0.925046297100680	,-0.921278554333542	,-0.917422293654362	,-0.913477885578796	,-0.909445709091850	,-0.905326151611455	,-0.901119608951253	,-0.896826485282557	,-0.892447193095527	,-0.887982153159530	,-0.883431794482714	,-0.878796554270789	,-0.874076877885018	,-0.869273218799429	,-0.864386038557238	,-0.859415806726511	,-0.854363000855042	,-0.849228106424471	,-0.844011616803638	,-0.838714033201178	,-0.833335864617367	,-0.827877627795215	,-0.822339847170813	,-0.816723054822951	,-0.811027790421990	,-0.805254601178011	,-0.799404041788239	,-0.793476674383745	,-0.787473068475438	,-0.781393800899343	,-0.775239455761179	,-0.769010624380237	,-0.762707905232565	,-0.756331903893465	,-0.749883232979311	,-0.743362512088683	,-0.736770367742838	,-0.730107433325514	,-0.723374349022072	,-0.716571761757984	,-0.709700325136678	,-0.702760699376740	,-0.695753551248474	,-0.688679554009844	,-0.681539387341780	,-0.674333737282879	,-0.667063296163485	,-0.659728762539171	,-0.652330841123620	,-0.644870242720914	,-0.637347684157243	,-0.629763888212025	,-0.622119583548463	,-0.614415504643537	,-0.606652391717429	,-0.598830990662405	,-0.590952052971147	,-0.583016335664548	,-0.575024601218978	,-0.566977617493023	,-0.558876157653708	,-0.550721000102209	,-0.542512928399065	,-0.534252731188889	,-0.525941202124599	,-0.517579139791156	,-0.509167347628841	,-0.500706633856054	,-0.492197811391662	,-0.483641697776893	,-0.475039115096782	,-0.466390889901188	,-0.457697853125374	,-0.448960840010173	,-0.440180690021732	,-0.431358246770859	,-0.422494357931968	,-0.413589875161628	,-0.404645654016739	,-0.395662553872328	,-0.386641437838978	,-0.377583172679897	,-0.368488628727644	,-0.359358679800500	,-0.350194203118512	,-0.340996079219210	,-0.331765191873003	,-0.322502427998261	,-0.313208677576105	,-0.303884833564891	,-0.294531791814416	,-0.285150450979842	,-0.275741712435352	,-0.266306480187545	,-0.256845660788576	,-0.247360163249056	,-0.237850898950711	,-0.228318781558813	,-0.218764726934398	,-0.209189653046265	,-0.199594479882777	,-0.189980129363470	,-0.180347525250468	,-0.170697593059730	,-0.161031259972125	,-0.151349454744346	,-0.141653107619675	,-0.131943150238601	,-0.122220515549308	,-0.112486137718038	,-0.102740952039331	,-0.0929858948461623	,-0.0832219034199790	,-0.0734499159006431	,-0.0636708711962937	,-0.0538857088931359	,-0.0440953691651631	,-0.0343007926838242	,-0.0245029205276421	,-0.0147026940917931	,-0.00490105499765685	,0.00490105499765685	,0.0147026940917931	,0.0245029205276421	,0.0343007926838242	,0.0440953691651631	,0.0538857088931359	,0.0636708711962937	,0.0734499159006431	,0.0832219034199790	,0.0929858948461623	,0.102740952039331	,0.112486137718038	,0.122220515549308	,0.131943150238601	,0.141653107619675	,0.151349454744346	,0.161031259972125	,0.170697593059730	,0.180347525250468	,0.189980129363470	,0.199594479882777	,0.209189653046265	,0.218764726934398	,0.228318781558813	,0.237850898950711	,0.247360163249056	,0.256845660788576	,0.266306480187545	,0.275741712435352	,0.285150450979842	,0.294531791814416	,0.303884833564891	,0.313208677576105	,0.322502427998261	,0.331765191873003	,0.340996079219210	,0.350194203118512	,0.359358679800500	,0.368488628727644	,0.377583172679897	,0.386641437838978	,0.395662553872328	,0.404645654016739	,0.413589875161628	,0.422494357931968	,0.431358246770859	,0.440180690021732	,0.448960840010173	,0.457697853125374	,0.466390889901188	,0.475039115096782	,0.483641697776893	,0.492197811391662	,0.500706633856054	,0.509167347628841	,0.517579139791156	,0.525941202124599	,0.534252731188889	,0.542512928399065	,0.550721000102209	,0.558876157653708	,0.566977617493023	,0.575024601218978	,0.583016335664548	,0.590952052971147	,0.598830990662405	,0.606652391717429	,0.614415504643537	,0.622119583548463	,0.629763888212025	,0.637347684157243	,0.644870242720914	,0.652330841123620	,0.659728762539171	,0.667063296163485	,0.674333737282879	,0.681539387341780	,0.688679554009844	,0.695753551248474	,0.702760699376740	,0.709700325136678	,0.716571761757984	,0.723374349022072	,0.730107433325514	,0.736770367742838	,0.743362512088683	,0.749883232979311	,0.756331903893465	,0.762707905232565	,0.769010624380237	,0.775239455761179	,0.781393800899343	,0.787473068475438	,0.793476674383745	,0.799404041788239	,0.805254601178011	,0.811027790421990	,0.816723054822951	,0.822339847170813	,0.827877627795215	,0.833335864617367	,0.838714033201178	,0.844011616803638	,0.849228106424471	,0.854363000855042	,0.859415806726511	,0.864386038557238	,0.869273218799429	,0.874076877885018	,0.878796554270789	,0.883431794482714	,0.887982153159530	,0.892447193095527	,0.896826485282557	,0.901119608951253	,0.905326151611455	,0.909445709091850	,0.913477885578796	,0.917422293654362	,0.921278554333542	,0.925046297100680	,0.928725159945059	,0.932314789395691	,0.935814840555278	,0.939224977133349	,0.942544871478572	,0.945774204610239	,0.948912666248912	,0.951959954846234	,0.954915777613911	,0.957779850551837	,0.960551898475386	,0.963231655041857	,0.965818862776063	,0.968313273095077	,0.970714646332118	,0.973022751759588	,0.975237367611243	,0.977358281103520	,0.979385288455988	,0.981318194910950	,0.983156814752184	,0.984900971322824	,0.986550497042385	,0.988105233422948	,0.989565031084503	,0.990929749769477	,0.992199258356498	,0.993373434873440	,0.994452166509907	,0.995435349629393	,0.996322889781677	,0.997114701716636	,0.997810709402343	,0.998410846055053	,0.998915054203942	,0.999323285871991	,0.999635503240465	,0.999851682193268	,0.999971849980221 };
		double GaussWeight[320] = { 7.22417022894800e-05	,0.000168158195616977	,0.000264200571979921	,0.000360229901104000	,0.000456227095866212	,0.000552181197609087	,0.000648082526890592	,0.000743921712961121	,0.000839689484512361	,0.000935376611415825	,0.00103097388562789	,0.00112647211435036	,0.00122186211760014	,0.00131713472754680	,0.00141228078862250	,0.00150729115799508	,0.00160215670622109	,0.00169686831799277	,0.00179141689293615	,0.00188579334643737	,0.00197998861048504	,0.00207399363452170	,0.00216779938630022	,0.00226139685274270	,0.00235477704080025	,0.00244793097831277	,0.00254084971486788	,0.00263352432265870	,0.00272594589733983	,0.00281810555888205	,0.00290999445242417	,0.00300160374912352	,0.00309292464700346	,0.00318394837179882	,0.00327466617779862	,0.00336506934868595	,0.00345514919837554	,0.00354489707184798	,0.00363430434598133	,0.00372336243037928	,0.00381206276819674	,0.00390039683696184	,0.00398835614939460	,0.00407593225422257	,0.00416311673699262	,0.00424990122087953	,0.00433627736749068	,0.00442223687766737	,0.00450777149228196	,0.00459287299303172	,0.00467753320322810	,0.00476174398858262	,0.00484549725798829	,0.00492878496429704	,0.00501159910509285	,0.00509393172346079	,0.00517577490875125	,0.00525712079734025	,0.00533796157338491	,0.00541828946957435	,0.00549809676787602	,0.00557737580027734	,0.00565611894952229	,0.00573431864984339	,0.00581196738768866	,0.00588905770244341	,0.00596558218714729	,0.00604153348920560	,0.00611690431109610	,0.00619168741106999	,0.00626587560384771	,0.00633946176130938	,0.00641243881317953	,0.00648479974770675	,0.00655653761233694	,0.00662764551438164	,0.00669811662168027	,0.00676794416325640	,0.00683712142996854	,0.00690564177515453	,0.00697349861527032	,0.00704068543052254	,0.00710719576549484	,0.00717302322976816	,0.00723816149853475	,0.00730260431320591	,0.00736634548201322	,0.00742937888060353	,0.00749169845262747	,0.00755329821032115	,0.00761417223508154	,0.00767431467803537	,0.00773371976060089	,0.00779238177504297	,0.00785029508502171	,0.00790745412613398	,0.00796385340644799	,0.00801948750703088	,0.00807435108246970	,0.00812843886138452	,0.00818174564693538	,0.00823426631732136	,0.00828599582627264	,0.00833692920353552	,0.00838706155535003	,0.00843638806491972	,0.00848490399287498	,0.00853260467772813	,0.00857948553632123	,0.00862554206426675	,0.00867076983638004	,0.00871516450710454	,0.00875872181092964	,0.00880143756280007	,0.00884330765851828	,0.00888432807513878	,0.00892449487135451	,0.00896380418787562	,0.00900225224780043	,0.00903983535697784	,0.00907654990436303	,0.00911239236236356	,0.00914735928717903	,0.00918144731913153	,0.00921465318298860	,0.00924697368827789	,0.00927840572959377	,0.00930894628689557	,0.00933859242579787	,0.00936734129785239	,0.00939519014082161	,0.00942213627894431	,0.00944817712319263	,0.00947331017152066	,0.00949753300910515	,0.00952084330857723	,0.00954323883024607	,0.00956471742231438	,0.00958527702108462	,0.00960491565115796	,0.00962363142562337	,0.00964142254623945	,0.00965828730360692	,0.00967422407733311	,0.00968923133618736	,0.00970330763824821	,0.00971645163104225	,0.00972866205167361	,0.00973993772694564	,0.00975027757347348	,0.00975968059778823	,0.00976814589643227	,0.00977567265604627	,0.00978226015344706	,0.00978790775569747	,0.00979261492016686	,0.00979638119458340	,0.00979920621707736	,0.00980108971621607	,0.00980203151103005	,0.00980203151103005	,0.00980108971621607	,0.00979920621707736	,0.00979638119458340	,0.00979261492016686	,0.00978790775569747	,0.00978226015344706	,0.00977567265604627	,0.00976814589643227	,0.00975968059778823	,0.00975027757347348	,0.00973993772694564	,0.00972866205167361	,0.00971645163104225	,0.00970330763824821	,0.00968923133618736	,0.00967422407733311	,0.00965828730360692	,0.00964142254623945	,0.00962363142562337	,0.00960491565115796	,0.00958527702108462	,0.00956471742231438	,0.00954323883024607	,0.00952084330857723	,0.00949753300910515	,0.00947331017152066	,0.00944817712319263	,0.00942213627894431	,0.00939519014082161	,0.00936734129785239	,0.00933859242579787	,0.00930894628689557	,0.00927840572959377	,0.00924697368827789	,0.00921465318298860	,0.00918144731913153	,0.00914735928717903	,0.00911239236236356	,0.00907654990436303	,0.00903983535697784	,0.00900225224780043	,0.00896380418787562	,0.00892449487135451	,0.00888432807513878	,0.00884330765851828	,0.00880143756280007	,0.00875872181092964	,0.00871516450710454	,0.00867076983638004	,0.00862554206426675	,0.00857948553632123	,0.00853260467772813	,0.00848490399287498	,0.00843638806491972	,0.00838706155535003	,0.00833692920353552	,0.00828599582627264	,0.00823426631732136	,0.00818174564693538	,0.00812843886138452	,0.00807435108246970	,0.00801948750703088	,0.00796385340644799	,0.00790745412613398	,0.00785029508502171	,0.00779238177504297	,0.00773371976060089	,0.00767431467803537	,0.00761417223508154	,0.00755329821032115	,0.00749169845262747	,0.00742937888060353	,0.00736634548201322	,0.00730260431320591	,0.00723816149853475	,0.00717302322976816	,0.00710719576549484	,0.00704068543052254	,0.00697349861527032	,0.00690564177515453	,0.00683712142996854	,0.00676794416325640	,0.00669811662168027	,0.00662764551438164	,0.00655653761233694	,0.00648479974770675	,0.00641243881317953	,0.00633946176130938	,0.00626587560384771	,0.00619168741106999	,0.00611690431109610	,0.00604153348920560	,0.00596558218714729	,0.00588905770244341	,0.00581196738768866	,0.00573431864984339	,0.00565611894952229	,0.00557737580027734	,0.00549809676787602	,0.00541828946957435	,0.00533796157338491	,0.00525712079734025	,0.00517577490875125	,0.00509393172346079	,0.00501159910509285	,0.00492878496429704	,0.00484549725798829	,0.00476174398858262	,0.00467753320322810	,0.00459287299303172	,0.00450777149228196	,0.00442223687766737	,0.00433627736749068	,0.00424990122087953	,0.00416311673699262	,0.00407593225422257	,0.00398835614939460	,0.00390039683696184	,0.00381206276819674	,0.00372336243037928	,0.00363430434598133	,0.00354489707184798	,0.00345514919837554	,0.00336506934868595	,0.00327466617779862	,0.00318394837179882	,0.00309292464700346	,0.00300160374912352	,0.00290999445242417	,0.00281810555888205	,0.00272594589733983	,0.00263352432265870	,0.00254084971486788	,0.00244793097831277	,0.00235477704080025	,0.00226139685274270	,0.00216779938630022	,0.00207399363452170	,0.00197998861048504	,0.00188579334643737	,0.00179141689293615	,0.00169686831799277	,0.00160215670622109	,0.00150729115799508	,0.00141228078862250	,0.00131713472754680	,0.00122186211760014	,0.00112647211435036	,0.00103097388562789	,0.000935376611415825	,0.000839689484512361	,0.000743921712961121	,0.000648082526890592	,0.000552181197609087	,0.000456227095866212	,0.000360229901104000	,0.000264200571979921	,0.000168158195616977	,7.22417022894800e-05 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

	//Gauss 480
	if (N == 480)
	{
		double Gauss[480] = { -0.999987475819599	,-0.999934011427531	,-0.999837827281408	,-0.999698906128460	,-0.999517251290280	,-0.999292869907203	,-0.999025771362748	,-0.998715966989650	,-0.998363469991356	,-0.997968295415516	,-0.997530460143254	,-0.997049982884035	,-0.996526884172755	,-0.995961186367790	,-0.995352913649438	,-0.994702092018559	,-0.994008749295255	,-0.993272915117555	,-0.992494620940076	,-0.991673900032619	,-0.990810787478714	,-0.989905320174096	,-0.988957536825114	,-0.987967477947057	,-0.986935185862424	,-0.985860704699096	,-0.984744080388457	,-0.983585360663418	,-0.982384595056378	,-0.981141834897104	,-0.979857133310535	,-0.978530545214509	,-0.977162127317415	,-0.975751938115769	,-0.974300037891710	,-0.972806488710426	,-0.971271354417499	,-0.969694700636173	,-0.968076594764554	,-0.966417105972724	,-0.964716305199783	,-0.962974265150823	,-0.961191060293812	,-0.959366766856416	,-0.957501462822737	,-0.955595227929984	,-0.953648143665058	,-0.951660293261073	,-0.949631761693796	,-0.947562635678016	,-0.945453003663836	,-0.943302955832892	,-0.941112584094499	,-0.938881982081718	,-0.936611245147361	,-0.934300470359909	,-0.931949756499361	,-0.929559204053018	,-0.927128915211184	,-0.924658993862793	,-0.922149545590978	,-0.919600677668547	,-0.917012499053407	,-0.914385120383897	,-0.911718653974067	,-0.909013213808870	,-0.906268915539293	,-0.903485876477414	,-0.900664215591384	,-0.897804053500344	,-0.894905512469268	,-0.891968716403735	,-0.888993790844636	,-0.885980862962803	,-0.882930061553577	,-0.879841517031298	,-0.876715361423734	,-0.873551728366435	,-0.870350753097020	,-0.867112572449398	,-0.863837324847917	,-0.860525150301447	,-0.857176190397395	,-0.853790588295653	,-0.850368488722478	,-0.846910037964304	,-0.843415383861492	,-0.839884675802006	,-0.836318064715029	,-0.832715703064511	,-0.829077744842655	,-0.825404345563325	,-0.821695662255410	,-0.817951853456101	,-0.814173079204121	,-0.810359501032879	,-0.806511281963569	,-0.802628586498198	,-0.798711580612555	,-0.794760431749117	,-0.790775308809888	,-0.786756382149184	,-0.782703823566346	,-0.778617806298400	,-0.774498505012646	,-0.770346095799197	,-0.766160756163451	,-0.761942665018497	,-0.757692002677476	,-0.753408950845866	,-0.749093692613717	,-0.744746412447825	,-0.740367296183845	,-0.735956531018349	,-0.731514305500821	,-0.727040809525601	,-0.722536234323762	,-0.718000772454940	,-0.713434617799101	,-0.708837965548252	,-0.704211012198100	,-0.699553955539649	,-0.694866994650745	,-0.690150329887569	,-0.685404162876070	,-0.680628696503347	,-0.675824134908973	,-0.670990683476273	,-0.666128548823543	,-0.661237938795215	,-0.656319062452975	,-0.651372130066825	,-0.646397353106095	,-0.641394944230402	,-0.636365117280562	,-0.631308087269445	,-0.626224070372786	,-0.621113283919945	,-0.615975946384614	,-0.610812277375481	,-0.605622497626838	,-0.600406828989149	,-0.595165494419566	,-0.589898717972396	,-0.584606724789523	,-0.579289741090789	,-0.573947994164316	,-0.568581712356796	,-0.563191125063727	,-0.557776462719606	,-0.552337956788081	,-0.546875839752057	,-0.541390345103754	,-0.535881707334728	,-0.530350161925849	,-0.524795945337231	,-0.519219294998127	,-0.513620449296778	,-0.507999647570223	,-0.502357130094067	,-0.496693138072211	,-0.491007913626540	,-0.485301699786574	,-0.479574740479079	,-0.473827280517639	,-0.468059565592191	,-0.462271842258521	,-0.456464357927728	,-0.450637360855646	,-0.444791100132228	,-0.438925825670904	,-0.433041788197894	,-0.427139239241492	,-0.421218431121312	,-0.415279616937503	,-0.409323050559930	,-0.403348986617321	,-0.397357680486381	,-0.391349388280878	,-0.385324366840693	,-0.379282873720842	,-0.373225167180464	,-0.367151506171782	,-0.361062150329036	,-0.354957359957379	,-0.348837396021754	,-0.342702520135738	,-0.336552994550356	,-0.330389082142872	,-0.324211046405554	,-0.318019151434404	,-0.311813661917877	,-0.305594843125558	,-0.299362960896829	,-0.293118281629499	,-0.286861072268423	,-0.280591600294084	,-0.274310133711164	,-0.268016941037084	,-0.261712291290527	,-0.255396453979938	,-0.249069699092004	,-0.242732297080112	,-0.236384518852785	,-0.230026635762108	,-0.223658919592121	,-0.217281642547206	,-0.210895077240449	,-0.204499496681985	,-0.198095174267332	,-0.191682383765699	,-0.185261399308285	,-0.178832495376561	,-0.172395946790538	,-0.165952028697014	,-0.159501016557820	,-0.153043186138035	,-0.146578813494207	,-0.140108174962546	,-0.133631547147115	,-0.127149206908003	,-0.120661431349492	,-0.114168497808210	,-0.107670683841277	,-0.101168267214439	,-0.0946615258901948	,-0.0881507380159124	,-0.0816361819119412	,-0.0751181360597126	,-0.0685968790898368	,-0.0620726897701908	,-0.0555458469940030	,-0.0490166297679303	,-0.0424853172001316	,-0.0359521884883366	,-0.0294175229079108	,-0.0228815997999171	,-0.0163446985591744	,-0.00980709862231460	,-0.00326907945583724	,0.00326907945583724	,0.00980709862231460	,0.0163446985591744	,0.0228815997999171	,0.0294175229079108	,0.0359521884883366	,0.0424853172001316	,0.0490166297679303	,0.0555458469940030	,0.0620726897701908	,0.0685968790898368	,0.0751181360597126	,0.0816361819119412	,0.0881507380159124	,0.0946615258901948	,0.101168267214439	,0.107670683841277	,0.114168497808210	,0.120661431349492	,0.127149206908003	,0.133631547147115	,0.140108174962546	,0.146578813494207	,0.153043186138035	,0.159501016557820	,0.165952028697014	,0.172395946790538	,0.178832495376561	,0.185261399308285	,0.191682383765699	,0.198095174267332	,0.204499496681985	,0.210895077240449	,0.217281642547206	,0.223658919592121	,0.230026635762108	,0.236384518852785	,0.242732297080112	,0.249069699092004	,0.255396453979938	,0.261712291290527	,0.268016941037084	,0.274310133711164	,0.280591600294084	,0.286861072268423	,0.293118281629499	,0.299362960896829	,0.305594843125558	,0.311813661917877	,0.318019151434404	,0.324211046405554	,0.330389082142872	,0.336552994550356	,0.342702520135738	,0.348837396021754	,0.354957359957379	,0.361062150329036	,0.367151506171782	,0.373225167180464	,0.379282873720842	,0.385324366840693	,0.391349388280878	,0.397357680486381	,0.403348986617321	,0.409323050559930	,0.415279616937503	,0.421218431121312	,0.427139239241492	,0.433041788197894	,0.438925825670904	,0.444791100132228	,0.450637360855646	,0.456464357927728	,0.462271842258521	,0.468059565592191	,0.473827280517639	,0.479574740479079	,0.485301699786574	,0.491007913626540	,0.496693138072211	,0.502357130094067	,0.507999647570223	,0.513620449296778	,0.519219294998127	,0.524795945337231	,0.530350161925849	,0.535881707334728	,0.541390345103754	,0.546875839752057	,0.552337956788081	,0.557776462719606	,0.563191125063727	,0.568581712356796	,0.573947994164316	,0.579289741090789	,0.584606724789523	,0.589898717972396	,0.595165494419566	,0.600406828989149	,0.605622497626838	,0.610812277375481	,0.615975946384614	,0.621113283919945	,0.626224070372786	,0.631308087269445	,0.636365117280562	,0.641394944230402	,0.646397353106095	,0.651372130066825	,0.656319062452975	,0.661237938795215	,0.666128548823543	,0.670990683476273	,0.675824134908973	,0.680628696503347	,0.685404162876070	,0.690150329887569	,0.694866994650745	,0.699553955539649	,0.704211012198100	,0.708837965548252	,0.713434617799101	,0.718000772454940	,0.722536234323762	,0.727040809525601	,0.731514305500821	,0.735956531018349	,0.740367296183845	,0.744746412447825	,0.749093692613717	,0.753408950845866	,0.757692002677476	,0.761942665018497	,0.766160756163451	,0.770346095799197	,0.774498505012646	,0.778617806298400	,0.782703823566346	,0.786756382149184	,0.790775308809888	,0.794760431749117	,0.798711580612555	,0.802628586498198	,0.806511281963569	,0.810359501032879	,0.814173079204121	,0.817951853456101	,0.821695662255410	,0.825404345563325	,0.829077744842655	,0.832715703064511	,0.836318064715029	,0.839884675802006	,0.843415383861492	,0.846910037964304	,0.850368488722478	,0.853790588295653	,0.857176190397395	,0.860525150301447	,0.863837324847917	,0.867112572449398	,0.870350753097020	,0.873551728366435	,0.876715361423734	,0.879841517031298	,0.882930061553577	,0.885980862962803	,0.888993790844636	,0.891968716403735	,0.894905512469268	,0.897804053500344	,0.900664215591384	,0.903485876477414	,0.906268915539293	,0.909013213808870	,0.911718653974067	,0.914385120383897	,0.917012499053407	,0.919600677668547	,0.922149545590978	,0.924658993862793	,0.927128915211184	,0.929559204053018	,0.931949756499361	,0.934300470359909	,0.936611245147361	,0.938881982081718	,0.941112584094499	,0.943302955832892	,0.945453003663836	,0.947562635678016	,0.949631761693796	,0.951660293261073	,0.953648143665058	,0.955595227929984	,0.957501462822737	,0.959366766856416	,0.961191060293812	,0.962974265150823	,0.964716305199783	,0.966417105972724	,0.968076594764554	,0.969694700636173	,0.971271354417499	,0.972806488710426	,0.974300037891710	,0.975751938115769	,0.977162127317415	,0.978530545214509	,0.979857133310535	,0.981141834897104	,0.982384595056378	,0.983585360663418	,0.984744080388457	,0.985860704699096	,0.986935185862424	,0.987967477947057	,0.988957536825114	,0.989905320174096	,0.990810787478714	,0.991673900032619	,0.992494620940076	,0.993272915117555	,0.994008749295255	,0.994702092018559	,0.995352913649438	,0.995961186367790	,0.996526884172755	,0.997049982884035	,0.997530460143254	,0.997968295415516	,0.998363469991356	,0.998715966989650	,0.999025771362748	,0.999292869907203	,0.999517251290280	,0.999698906128460	,0.999837827281408	,0.999934011427531	,0.999987475819599 };
		double GaussWeight[480] = { 3.21410242465769e-05	,7.48168530395741e-05	,0.000117552677581628	,0.000160288965945794	,0.000203019504362265	,0.000245741694281466	,0.000288453504781483	,0.000331153040461037	,0.000373838448026151	,0.000416507890053588	,0.000459159536180179	,0.000501791559744907	,0.000544402136392781	,0.000586989443465939	,0.000629551659738542	,0.000672086965313827	,0.000714593541601561	,0.000757069571337964	,0.000799513238628749	,0.000841922729005419	,0.000884296229489298	,0.000926631928660088	,0.000968928016727617	,0.00101118268560508	,0.00105339412898368	,0.00109556054240786	,0.00113768012335108	,0.00117975107129187	,0.00122177158779008	,0.00126373987656320	,0.00130565414356283	,0.00134751259705092	,0.00138931344767620	,0.00143105490855053	,0.00147273519532516	,0.00151435252626673	,0.00155590512233356	,0.00159739120725153	,0.00163880900758991	,0.00168015675283727	,0.00172143267547693	,0.00176263501106276	,0.00180376199829422	,0.00184481187909201	,0.00188578289867290	,0.00192667330562497	,0.00196748135198231	,0.00200820529329971	,0.00204884338872744	,0.00208939390108541	,0.00212985509693755	,0.00217022524666589	,0.00221050262454449	,0.00225068550881319	,0.00229077218175123	,0.00233076092975072	,0.00237065004338973	,0.00241043781750555	,0.00245012255126747	,0.00248970254824958	,0.00252917611650311	,0.00256854156862893	,0.00260779722184964	,0.00264694139808145	,0.00268597242400596	,0.00272488863114167	,0.00276368835591528	,0.00280236993973290	,0.00284093172905083	,0.00287937207544633	,0.00291768933568797	,0.00295588187180608	,0.00299394805116259	,0.00303188624652087	,0.00306969483611526	,0.00310737220372051	,0.00314491673872071	,0.00318232683617834	,0.00321960089690262	,0.00325673732751811	,0.00329373454053270	,0.00333059095440545	,0.00336730499361438	,0.00340387508872344	,0.00344029967645008	,0.00347657719973169	,0.00351270610779240	,0.00354868485620917	,0.00358451190697808	,0.00362018572857976	,0.00365570479604511	,0.00369106759102034	,0.00372627260183193	,0.00376131832355135	,0.00379620325805913	,0.00383092591410916	,0.00386548480739240	,0.00389987846060013	,0.00393410540348728	,0.00396816417293536	,0.00400205331301472	,0.00403577137504707	,0.00406931691766719	,0.00410268850688477	,0.00413588471614551	,0.00416890412639214	,0.00420174532612521	,0.00423440691146326	,0.00426688748620298	,0.00429918566187873	,0.00433130005782208	,0.00436322930122066	,0.00439497202717700	,0.00442652687876673	,0.00445789250709668	,0.00448906757136249	,0.00452005073890601	,0.00455084068527214	,0.00458143609426556	,0.00461183565800691	,0.00464203807698881	,0.00467204206013127	,0.00470184632483692	,0.00473144959704599	,0.00476085061129060	,0.00479004811074880	,0.00481904084729853	,0.00484782758157081	,0.00487640708300278	,0.00490477812989015	,0.00493293950943977	,0.00496089001782108	,0.00498862846021786	,0.00501615365087913	,0.00504346441316988	,0.00507055957962144	,0.00509743799198147	,0.00512409850126302	,0.00515053996779422	,0.00517676126126672	,0.00520276126078395	,0.00522853885490910	,0.00525409294171282	,0.00527942242881994	,0.00530452623345655	,0.00532940328249602	,0.00535405251250497	,0.00537847286978878	,0.00540266331043661	,0.00542662280036587	,0.00545035031536665	,0.00547384484114543	,0.00549710537336828	,0.00552013091770410	,0.00554292048986677	,0.00556547311565763	,0.00558778783100664	,0.00560986368201403	,0.00563169972499089	,0.00565329502649939	,0.00567464866339300	,0.00569575972285563	,0.00571662730244081	,0.00573725051011031	,0.00575762846427208	,0.00577776029381818	,0.00579764513816187	,0.00581728214727436	,0.00583667048172127	,0.00585580931269850	,0.00587469782206749	,0.00589333520239045	,0.00591172065696467	,0.00592985339985670	,0.00594773265593587	,0.00596535766090744	,0.00598272766134534	,0.00599984191472425	,0.00601669968945148	,0.00603330026489816	,0.00604964293143001	,0.00606572699043789	,0.00608155175436729	,0.00609711654674802	,0.00611242070222309	,0.00612746356657709	,0.00614224449676417	,0.00615676286093544	,0.00617101803846616	,0.00618500941998220	,0.00619873640738592	,0.00621219841388208	,0.00622539486400265	,0.00623832519363147	,0.00625098885002829	,0.00626338529185271	,0.00627551398918692	,0.00628737442355866	,0.00629896608796316	,0.00631028848688493	,0.00632134113631903	,0.00633212356379152	,0.00634263530837985	,0.00635287592073249	,0.00636284496308820	,0.00637254200929463	,0.00638196664482665	,0.00639111846680402	,0.00639999708400858	,0.00640860211690103	,0.00641693319763719	,0.00642498997008357	,0.00643277208983275	,0.00644027922421806	,0.00644751105232773	,0.00645446726501868	,0.00646114756492972	,0.00646755166649441	,0.00647367929595283	,0.00647953019136380	,0.00648510410261572	,0.00649040079143745	,0.00649542003140844	,0.00650016160796839	,0.00650462531842645	,0.00650881097196977	,0.00651271838967192	,0.00651634740450028	,0.00651969786132324	,0.00652276961691705	,0.00652556253997146	,0.00652807651109590	,0.00653031142282418	,0.00653226717961933	,0.00653394369787742	,0.00653534090593135	,0.00653645874405399	,0.00653729716446038	,0.00653785613131014	,0.00653813562070871	,0.00653813562070871	,0.00653785613131014	,0.00653729716446038	,0.00653645874405399	,0.00653534090593135	,0.00653394369787742	,0.00653226717961933	,0.00653031142282418	,0.00652807651109590	,0.00652556253997146	,0.00652276961691705	,0.00651969786132324	,0.00651634740450028	,0.00651271838967192	,0.00650881097196977	,0.00650462531842645	,0.00650016160796839	,0.00649542003140844	,0.00649040079143745	,0.00648510410261572	,0.00647953019136380	,0.00647367929595283	,0.00646755166649441	,0.00646114756492972	,0.00645446726501868	,0.00644751105232773	,0.00644027922421806	,0.00643277208983275	,0.00642498997008357	,0.00641693319763719	,0.00640860211690103	,0.00639999708400858	,0.00639111846680402	,0.00638196664482665	,0.00637254200929463	,0.00636284496308820	,0.00635287592073249	,0.00634263530837985	,0.00633212356379152	,0.00632134113631903	,0.00631028848688493	,0.00629896608796316	,0.00628737442355866	,0.00627551398918692	,0.00626338529185271	,0.00625098885002829	,0.00623832519363147	,0.00622539486400265	,0.00621219841388208	,0.00619873640738592	,0.00618500941998220	,0.00617101803846616	,0.00615676286093544	,0.00614224449676417	,0.00612746356657709	,0.00611242070222309	,0.00609711654674802	,0.00608155175436729	,0.00606572699043789	,0.00604964293143001	,0.00603330026489816	,0.00601669968945148	,0.00599984191472425	,0.00598272766134534	,0.00596535766090744	,0.00594773265593587	,0.00592985339985670	,0.00591172065696467	,0.00589333520239045	,0.00587469782206749	,0.00585580931269850	,0.00583667048172127	,0.00581728214727436	,0.00579764513816187	,0.00577776029381818	,0.00575762846427208	,0.00573725051011031	,0.00571662730244081	,0.00569575972285563	,0.00567464866339300	,0.00565329502649939	,0.00563169972499089	,0.00560986368201403	,0.00558778783100664	,0.00556547311565763	,0.00554292048986677	,0.00552013091770410	,0.00549710537336828	,0.00547384484114543	,0.00545035031536665	,0.00542662280036587	,0.00540266331043661	,0.00537847286978878	,0.00535405251250497	,0.00532940328249602	,0.00530452623345655	,0.00527942242881994	,0.00525409294171282	,0.00522853885490910	,0.00520276126078395	,0.00517676126126672	,0.00515053996779422	,0.00512409850126302	,0.00509743799198147	,0.00507055957962144	,0.00504346441316988	,0.00501615365087913	,0.00498862846021786	,0.00496089001782108	,0.00493293950943977	,0.00490477812989015	,0.00487640708300278	,0.00484782758157081	,0.00481904084729853	,0.00479004811074880	,0.00476085061129060	,0.00473144959704599	,0.00470184632483692	,0.00467204206013127	,0.00464203807698881	,0.00461183565800691	,0.00458143609426556	,0.00455084068527214	,0.00452005073890601	,0.00448906757136249	,0.00445789250709668	,0.00442652687876673	,0.00439497202717700	,0.00436322930122066	,0.00433130005782208	,0.00429918566187873	,0.00426688748620298	,0.00423440691146326	,0.00420174532612521	,0.00416890412639214	,0.00413588471614551	,0.00410268850688477	,0.00406931691766719	,0.00403577137504707	,0.00400205331301472	,0.00396816417293536	,0.00393410540348728	,0.00389987846060013	,0.00386548480739240	,0.00383092591410916	,0.00379620325805913	,0.00376131832355135	,0.00372627260183193	,0.00369106759102034	,0.00365570479604511	,0.00362018572857976	,0.00358451190697808	,0.00354868485620917	,0.00351270610779240	,0.00347657719973169	,0.00344029967645008	,0.00340387508872344	,0.00336730499361438	,0.00333059095440545	,0.00329373454053270	,0.00325673732751811	,0.00321960089690262	,0.00318232683617834	,0.00314491673872071	,0.00310737220372051	,0.00306969483611526	,0.00303188624652087	,0.00299394805116259	,0.00295588187180608	,0.00291768933568797	,0.00287937207544633	,0.00284093172905083	,0.00280236993973290	,0.00276368835591528	,0.00272488863114167	,0.00268597242400596	,0.00264694139808145	,0.00260779722184964	,0.00256854156862893	,0.00252917611650311	,0.00248970254824958	,0.00245012255126747	,0.00241043781750555	,0.00237065004338973	,0.00233076092975072	,0.00229077218175123	,0.00225068550881319	,0.00221050262454449	,0.00217022524666589	,0.00212985509693755	,0.00208939390108541	,0.00204884338872744	,0.00200820529329971	,0.00196748135198231	,0.00192667330562497	,0.00188578289867290	,0.00184481187909201	,0.00180376199829422	,0.00176263501106276	,0.00172143267547693	,0.00168015675283727	,0.00163880900758991	,0.00159739120725153	,0.00155590512233356	,0.00151435252626673	,0.00147273519532516	,0.00143105490855053	,0.00138931344767620	,0.00134751259705092	,0.00130565414356283	,0.00126373987656320	,0.00122177158779008	,0.00117975107129187	,0.00113768012335108	,0.00109556054240786	,0.00105339412898368	,0.00101118268560508	,0.000968928016727617	,0.000926631928660088	,0.000884296229489298	,0.000841922729005419	,0.000799513238628749	,0.000757069571337964	,0.000714593541601561	,0.000672086965313827	,0.000629551659738542	,0.000586989443465939	,0.000544402136392781	,0.000501791559744907	,0.000459159536180179	,0.000416507890053588	,0.000373838448026151	,0.000331153040461037	,0.000288453504781483	,0.000245741694281466	,0.000203019504362265	,0.000160288965945794	,0.000117552677581628	,7.48168530395741e-05	,3.21410242465769e-05 };
		memcpy(Point, Gauss, sizeof(Gauss));
		memcpy(Weight, GaussWeight, sizeof(GaussWeight));
	}

}

void HSMA2D::ConstructLeftTerm(double** LeftTermReal, double** LeftTermImag, double*** IntegralTop, double*** IntegralDown, int Nw, int NJKBound, int p, double EU, double mu, double Lx, double Ly, double Lz, int S)
{

	double Q[S * Nw][S * Nw][p * p], QZD[S * Nw][S * Nw][p * p];
	double DownQ[S * Nw][S * Nw][p * p], DownQZD[S * Nw][S * Nw][p * p];


#pragma omp parallel 
	{
		double QM[p * p], QZDM[p * p];
		int index_r, index_dr;
#pragma omp for schedule(static) private(QM,QZDM)
		for (int m = 0; m < S * Nw; m++)
			for (int n = 0; n < S * Nw; n++)
			{
				CalculateMultipoleExpansion(QM, p, IntegralTop[m][n][0], IntegralTop[m][n][1], IntegralTop[m][n][2]);
				CalculateZDerivativeMultipoleExpansion(QZDM, p, IntegralTop[m][n][0], IntegralTop[m][n][1], IntegralTop[m][n][2]);
				for (int i = 0; i < p * p; i++)
				{
					Q[m][n][i] = QM[i];
					QZD[m][n][i] = QZDM[i];
				}

				CalculateMultipoleExpansion(QM, p, IntegralDown[m][n][0], IntegralDown[m][n][1], IntegralDown[m][n][2]);
				CalculateZDerivativeMultipoleExpansion(QZDM, p, IntegralDown[m][n][0], IntegralDown[m][n][1], IntegralDown[m][n][2]);
				for (int i = 0; i < p * p; i++)
				{
					DownQ[m][n][i] = QM[i];
					DownQZD[m][n][i] = QZDM[i];
				}
			}

		double Vec1, Vec2, Vec3, Vec4, I1, I2, I3, I4, EQ1, EQ2, EQ3, EQ4;
		int Up, Down;

#pragma omp for schedule(static)
		for (int j = -NJKBound; j <= NJKBound; j++)
			for (int k = -NJKBound; k <= NJKBound; k++)
			{
				index_r = (j + NJKBound) * (2 * NJKBound + 1) + k + NJKBound;
				index_dr = (2 * NJKBound + 1) * (2 * NJKBound + 1) + (j + NJKBound) * (2 * NJKBound + 1) + k + NJKBound;
				for (int m = 0; m < S * Nw; m++)
					for (int n = 0; n < S * Nw; n++)
					{
						I1 = -mu * (j * IntegralTop[m][n][0] + k * IntegralTop[m][n][1]);
						for (int i = 1; i < p * p; i++)
						{
							Up = floor((p * p - 1) / 2);
							Down = ceil((p * p - 1) / 2.0);

							Vec1 = (1.0 / EU) * QZD[m][n][i] + mu * sqrt(j * j + k * k) * Q[m][n][i];

							if (i % 2 == 0)
							{
								EQ1 = Vec1 * cos(I1);
								LeftTermReal[index_r][i / 2 - 1] -= EQ1 * IntegralTop[m][n][3];
							}

							if (i % 2 == 1)
							{
								EQ1 = Vec1 * sin(I1);
								LeftTermImag[index_r][(i - 1) / 2] -= EQ1 * IntegralTop[m][n][3];
							}

							Vec1 = (1.0 / EU) * DownQZD[m][n][i] - mu * sqrt(j * j + k * k) * DownQ[m][n][i];

							if (i % 2 == 0)
							{
								EQ1 = Vec1 * cos(I1);
								LeftTermReal[index_dr][i / 2 - 1] -= EQ1 * IntegralTop[m][n][3];
							}

							if (i % 2 == 1)
							{
								EQ1 = Vec1 * sin(I1);
								LeftTermImag[index_dr][(i - 1) / 2] -= EQ1 * IntegralTop[m][n][3];
							}
						}
					}
			}

	}

}

double HSMA2D::fac(double t)//calculate factorial
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

void HSMA2D::CalculateRDMultipoleExpansion(double* Q, int p, double x, double y, double z)
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

void HSMA2D::CalculateLocalRDMultipoleExpansion(double* Q, int p, double x, double y, double z, double Rs)
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

void HSMA2D::CalculateMultipoleExpansion(double* Q, int p, double x, double y, double z)
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

void HSMA2D::CalculateZDerivativeMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
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

void HSMA2D::CalculateXDMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
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

void HSMA2D::CalculateYDMultipoleExpansion(double* Q, int p, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
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

void HSMA2D::SetFibonacci(double **Fibonacci, double F, double Fp, int Np, double Rs, double PI)
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

void HSMA2D::AdjustParticle_Double(double Particle[][3], int N, double Lx, double Ly, double Lz)//Image Correction -Lx/2~Lx/2
{
	for (int i = 0; i < N; i++)
	{
		if (Particle[i][0] > Lx / 2)
		{
			Particle[i][0] = Particle[i][0] - ceil((Particle[i][0] - Lx / 2.0) / (Lx + 0.00)) * Lx;
		}
		if (Particle[i][0] < -Lx / 2)
		{
			Particle[i][0] = Particle[i][0] + ceil(fabs((Particle[i][0] + Lx / 2.0) / Lx)) * Lx;
		}
		if (Particle[i][1] > Ly / 2)
		{
			Particle[i][1] = Particle[i][1] - ceil((Particle[i][1] - Ly / 2.0) / (Ly + 0.00)) * Ly;
		}
		if (Particle[i][1] < -Ly / 2)
		{
			Particle[i][1] = Particle[i][1] + ceil(fabs((Particle[i][1] + Ly / 2.0) / Ly)) * Ly;
		}
	}
}

double HSMA2D::SetImageCharge(double ImageCharge[][5], int* ImageNumber, int TotalNumber, double Source[][3], double* Q, int NSource, double Rs, double Lx, double Ly, double Lz, int lx, int ly, int lz)
{
	int total = 0;
	int number = 0;
	for (int i = -lx - 1; i <= lx + 1; i++)
	{
		double CX, CY, CZ;
		for (int j = -ly - 1; j <= ly + 1; j++)
			for (int k = -lz; k <= lz; k++)
			{
				for (int m = 0; m < NSource; m++)
				{
					CX = Source[m][0] + i * Lx;
					CY = Source[m][1] + j * Ly;
					CZ = pow(-1, k) * Source[m][2] + k * Lz;
					if (CX * CX + CY * CY + CZ * CZ <= Rs * Rs)
					{					   
						ImageCharge[number][0] = CX;
						ImageCharge[number][1] = CY;
						ImageCharge[number][2] = CZ;
						ImageCharge[number][3] = Q[m];
						ImageCharge[number][4] = k;
						number++;
					}
				}
			}
	}

	*ImageNumber = number;
	return 1.0;
}

void HSMA2D::CalculateNearFieldAndZD(double** Top, double** TopZD, double** Down, double** DownZD, double ImageCharge[][5], int ImageNumber, double*** IntegralTop, double*** IntegralDown, int Nw, double Gamma, int IF_FMM_RightTerm, int S, double* AR, double Lx, double Ly, double Lz, double tolerance)
{
	if (IF_FMM_RightTerm)//Using FMM to calculate pairwise sum
	{
		double eps = tolerance;

		/*            Set FMM parameters           */
		int ns = ImageNumber;
		int nt = (S * Nw) * (S * Nw) * 2;
		double* source = (double*)malloc(3 * ns * sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double* charge = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			source[3 * i] = ImageCharge[i][0];
			source[3 * i + 1] = ImageCharge[i][1];
			source[3 * i + 2] = ImageCharge[i][2];
			charge[i] = ImageCharge[i][3] * pow(Gamma, fabs(ImageCharge[i][4]) + 0.00);
		}

		for (int i = 0; i < S * Nw; i++)
			for (int j = 0; j < S * Nw; j++)
			{
				target[3 * (i * (S * Nw) + j)] = IntegralTop[i][j][0];
				target[3 * (i * (S * Nw) + j) + 1] = IntegralTop[i][j][1];
				target[3 * (i * (S * Nw) + j) + 2] = IntegralTop[i][j][2];

				target[3 * ((S * Nw) * (S * Nw) + i * (S * Nw) + j)] = IntegralDown[i][j][0];
				target[3 * ((S * Nw) * (S * Nw) + i * (S * Nw) + j) + 1] = IntegralDown[i][j][1];
				target[3 * ((S * Nw) * (S * Nw) + i * (S * Nw) + j) + 2] = IntegralDown[i][j][2];
			}

		double* pottarg = (double*)malloc(nt * sizeof(double));
		double* gradtarg = (double*)malloc(3 * nt * sizeof(double));

		int ier = 0;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg, &ier);

		for (int i = 0; i < S * Nw; i++)
			for (int j = 0; j < S * Nw; j++)
			{
				Top[i][j] = pottarg[i * (S * Nw) + j];
				TopZD[i][j] = gradtarg[3 * (i * (S * Nw) + j) + 2];
				Down[i][j] = pottarg[(S * Nw) * (S * Nw) + i * (S * Nw) + j];
				DownZD[i][j] = gradtarg[3 * ((S * Nw) * (S * Nw) + i * (S * Nw) + j) + 2];
			}
		 
		free(source); free(target); free(charge); free(pottarg); free(gradtarg);
	}
	else
	{
		double Paramet1[ImageNumber];
		for (int i = 0; i < ImageNumber; i++)
		{
			Paramet1[i] = ImageCharge[i][3] * pow(Gamma, fabs(ImageCharge[i][4]) + 0.00);
		}
		
        #pragma omp parallel
		{
			int id = omp_get_thread_num();
			int size = omp_get_num_threads();
			int min_atom=id*floor(S*Nw*S*Nw/size)+1, max_atom=(id+1)*floor(S*Nw*S*Nw/size);
			if (id == size - 1)max_atom = S * Nw * S * Nw - 1;
			if (id == 0)min_atom = 0;

			int float_double;
			if (tolerance > 0.000001) { 
				float_double = 1; 
				float IntegralTop_X[int(ceil( (max_atom - min_atom + 1)/16.0 ))*16], IntegralTop_Y[int(ceil((max_atom - min_atom + 1) / 16.0)) * 16], IntegralTop_Z[int(ceil((max_atom - min_atom + 1) / 16.0)) * 16];
				for (int i = min_atom; i <= max_atom; i++)
				{
					int ix = floor(i / (S * Nw));
					int iy = i - ix * (S * Nw);
					IntegralTop_X[i - min_atom] = IntegralTop[ix][iy][0];
					IntegralTop_Y[i - min_atom] = IntegralTop[ix][iy][1];
					IntegralTop_Z[i - min_atom] = IntegralTop[ix][iy][2];
				}

				for (int i = max_atom - min_atom + 1; i<int(ceil((max_atom - min_atom + 1) / 16.0)) * 16;i++)
				{
					IntegralTop_X[i] = 0.00;
					IntegralTop_Y[i] = 0.00;
					IntegralTop_Z[i] = 0.00;
				}

				for (int i = min_atom; i <= max_atom; i = i + 16)
				{
					__m512 pottarg, fldtarg, pottarg2, fldtarg2, X0, Y0, X1, Y1, dx, dy, dz, dz1, delta, delta1, Para, midterm, midterm1;
					float F1[16], F2[16], F3[16], F4[16];
					pottarg = fldtarg = pottarg2 = fldtarg2 = _mm512_setzero_ps();
					X0 = _mm512_load_ps(&IntegralTop_X[i - min_atom]);
					Y0 = _mm512_load_ps(&IntegralTop_Y[i - min_atom]);

					for (int j = 0; j < ImageNumber; j++)
					{
						Para = _mm512_set1_ps(Paramet1[j]);
						X1 = _mm512_set1_ps(ImageCharge[j][0]);
						Y1 = _mm512_set1_ps(ImageCharge[j][1]);
						dz = _mm512_set1_ps(Lz / 2.0 - ImageCharge[j][2]);
						dx = X0 - X1;
						dy = Y0 - Y1;
						delta = _mm512_invsqrt_ps(dx * dx + dy * dy + dz * dz);
						midterm = Para * delta;
						fldtarg -= dz * midterm * delta * delta;
						pottarg += midterm;
						dz1 = _mm512_set1_ps(-Lz / 2.0 - ImageCharge[j][2]);
						delta1= _mm512_invsqrt_ps(dx * dx + dy * dy + dz1 * dz1);
						midterm1 = Para * delta1;
						pottarg2 += midterm1;
						fldtarg2 -= dz1 * midterm1 * delta1 * delta1;
					}

					_mm512_store_ps(&F1[0], pottarg);
					_mm512_store_ps(&F2[0], pottarg2);
					_mm512_store_ps(&F4[0], fldtarg2);
					_mm512_store_ps(&F3[0], fldtarg);

					for (int j = i; (j < i + 16) && (j <= max_atom); j++)
					{
						int ix = floor(j / (S * Nw));
						int iy = j - ix * (S * Nw);

						Top[ix][iy] = F1[j - i];
						TopZD[ix][iy] = F3[j - i];
						Down[ix][iy] = F2[j - i];
						DownZD[ix][iy] = F4[j - i];
					}
				}
			}
			else { 
				float_double = 2; 
				double IntegralTop_X[int(ceil((max_atom - min_atom + 1) / 8.0)) * 8], IntegralTop_Y[int(ceil((max_atom - min_atom + 1) / 8.0)) * 8], IntegralTop_Z[int(ceil((max_atom - min_atom + 1) / 8.0)) * 8];
				for (int i = min_atom; i <= max_atom; i++)
				{
					int ix = floor(i / (S * Nw));
					int iy = i - ix * (S * Nw);
					IntegralTop_X[i - min_atom] = IntegralTop[ix][iy][0];
					IntegralTop_Y[i - min_atom] = IntegralTop[ix][iy][1];
					IntegralTop_Z[i - min_atom] = IntegralTop[ix][iy][2];
				}

				for (int i = max_atom - min_atom + 1; i<int(ceil((max_atom - min_atom + 1) / 8.0)) * 8; i++)
				{
					IntegralTop_X[i] = 0.00;
					IntegralTop_Y[i] = 0.00;
					IntegralTop_Z[i] = 0.00;
				}

				for (int i = min_atom; i <= max_atom; i = i + 8)
				{
					__m512d pottarg, fldtarg, pottarg2, fldtarg2, X0, Y0, X1, Y1, dx, dy, dz, dz1, delta, delta1, Para, midterm, midterm1;
					double F1[8], F2[8], F3[8], F4[8];
					pottarg = fldtarg = pottarg2 = fldtarg2 = _mm512_setzero_pd();
					X0 = _mm512_load_pd(&IntegralTop_X[i - min_atom]);
					Y0 = _mm512_load_pd(&IntegralTop_Y[i - min_atom]);

					for (int j = 0; j < ImageNumber; j++)
					{
						Para = _mm512_set1_pd(Paramet1[j]);
						X1 = _mm512_set1_pd(ImageCharge[j][0]);
						Y1 = _mm512_set1_pd(ImageCharge[j][1]);
						dz = _mm512_set1_pd(Lz / 2.0 - ImageCharge[j][2]);
						dx = X0 - X1;
						dy = Y0 - Y1;
						delta = _mm512_invsqrt_pd(dx * dx + dy * dy + dz * dz);
						midterm = Para * delta;
						pottarg += midterm;
						fldtarg -= dz * midterm * delta * delta;
						dz1 = _mm512_set1_pd(-Lz / 2.0 - ImageCharge[j][2]);
						delta1 = _mm512_invsqrt_pd(dx * dx + dy * dy + dz1 * dz1);
						midterm1 = Para * delta1;
						pottarg2 += midterm1;
						fldtarg2 -= dz1 * midterm1 * delta1 * delta1;
					}

					_mm512_store_pd(&F1[0], pottarg);
					_mm512_store_pd(&F2[0], pottarg2);
					_mm512_store_pd(&F4[0], fldtarg2);
					_mm512_store_pd(&F3[0], fldtarg);

					for (int j = i; (j < i + 8) && (j < max_atom); j++)
					{
						int ix = floor(j / (S * Nw));
						int iy = j - ix * (S * Nw);

						Top[ix][iy] = F1[j - i];
						TopZD[ix][iy] = F3[j - i];
						Down[ix][iy] = F2[j - i];
						DownZD[ix][iy] = F4[j - i];
					}
				}
			}
		}
		
	}
}

void HSMA2D::ConstructRightTerm(double* RightTermReal, double* RightTermImag, double** TopNear, double** TopZDNear, double** DownNear, double** DownZDNear, double*** IntegralTop, double*** IntegralDown, int Nw, int NJKBound, double EU, double mu, double Lx, double Ly, double Lz, int S, double tolerance)
{

	for (int i = 0; i < 2 * (2 * NJKBound + 1) * (2 * NJKBound + 1); i++)
	{
		RightTermReal[i] = 0.00;
		RightTermImag[i] = 0.00;
	}

     #pragma omp parallel
	{
		int id = omp_get_thread_num();
		int size = omp_get_num_threads();
		int min_index = id * floor((2 * NJKBound + 1) * (2 * NJKBound + 1) / (size+0.00)) + 1, max_index = (id + 1) * floor((2 * NJKBound + 1) * (2 * NJKBound + 1) / (size+0.00));
		if (id == size - 1)max_index = (2 * NJKBound + 1) * (2 * NJKBound + 1) - 1;
		if (id == 0)min_index = 0;

		if (tolerance > 0.000001)
		{
			float IntegralTop_X[int(ceil(S * Nw * S * Nw / 16.0)) * 16], IntegralTop_Y[int(ceil(S * Nw * S * Nw / 16.0)) * 16], IntegralTop_W[int(ceil(S * Nw * S * Nw / 16.0)) * 16];
			float IntegralDown_X[int(ceil(S * Nw * S * Nw / 16.0)) * 16], IntegralDown_Y[int(ceil(S * Nw * S * Nw / 16.0)) * 16], IntegralDown_W[int(ceil(S * Nw * S * Nw / 16.0)) * 16];
			float TopNear_F[int(ceil(S * Nw * S * Nw / 16.0)) * 16], TopZDNear_F[int(ceil(S * Nw * S * Nw / 16.0)) * 16], DownNear_F[int(ceil(S * Nw * S * Nw / 16.0)) * 16], DownZDNear_F[int(ceil(S * Nw * S * Nw / 16.0)) * 16];
			for (int i = 0; i < S * Nw * S * Nw; i++)
			{
				int ix = floor(i / (S * Nw));
				int iy = i - ix * (S * Nw);
				IntegralTop_X[i] = IntegralTop[ix][iy][0];
				IntegralTop_Y[i] = IntegralTop[ix][iy][1];
				IntegralTop_W[i] = IntegralTop[ix][iy][3];
				IntegralDown_X[i] = IntegralDown[ix][iy][0];
				IntegralDown_Y[i] = IntegralDown[ix][iy][1];
				IntegralDown_W[i] = IntegralDown[ix][iy][3];
				TopNear_F[i] = TopNear[ix][iy];
				TopZDNear_F[i] = TopZDNear[ix][iy];
				DownNear_F[i] = DownNear[ix][iy];
				DownZDNear_F[i] = DownZDNear[ix][iy];
			}

			for (int i = S * Nw * S * Nw; i < int(ceil(S * Nw * S * Nw / 16.0)) * 16; i++)
			{
				IntegralTop_X[i] = 0.00;
				IntegralTop_Y[i] = 0.00;
				IntegralDown_X[i] = 0.00;
				IntegralDown_Y[i] = 0.00;
				TopNear_F[i] = 0.00;
				TopZDNear_F[i] = 0.00;
				DownNear_F[i] = 0.00;
				DownZDNear_F[i] = 0.00;
				IntegralTop_W[i] = 0.00;
				IntegralDown_W[i] = 0.00;
			}

			float inveu = (1.0 / EU);

			for (int i = min_index; i <= max_index; i++)
			{
				int j = floor(i / ((2 * NJKBound + 1)));
				int k = i - j * ((2 * NJKBound + 1));
				j = j - NJKBound;
				k = k - NJKBound;

				__m512 X0, Y0, W0, topnear, topzdnear, downnear, downzdnear, I1, E1, E2, I2, Sin1, Sin2, Cos1, Cos2, EQ1, EQ2, EQ3, EQ4;
				__m512 MU, J, K, INVEU;

				MU = _mm512_set1_ps(mu);
				J = _mm512_set1_ps(j);
				K = _mm512_set1_ps(k);
				INVEU = _mm512_set1_ps(inveu);

				float res1 = 0.00, res2 = 0.00, res3 = 0.00, res4 = 0.00;

				for (int m = 0; m < S * Nw * S * Nw; m += 16)
				{
					X0 = _mm512_load_ps(&IntegralTop_X[m]);
					Y0 = _mm512_load_ps(&IntegralTop_Y[m]);
					W0 = _mm512_load_ps(&IntegralTop_W[m]);
					topnear = _mm512_load_ps(&TopNear_F[m]);
					topzdnear = _mm512_load_ps(&TopZDNear_F[m]);
					downnear = _mm512_load_ps(&DownNear_F[m]);
					downzdnear = _mm512_load_ps(&DownZDNear_F[m]);

					I1 = -MU * (J * X0 + K * Y0);
					E1 = INVEU * topzdnear + MU * _mm512_sqrt_ps(J * J + K * K) * topnear;
					Sin1 = _mm512_sincos_ps(&Cos1, I1);
					EQ1 = E1 * Cos1;
					EQ2 = E1 * Sin1;
					res1 += _mm512_reduce_add_ps(EQ1 * W0);
					res2 += _mm512_reduce_add_ps(EQ2 * W0);

					E2 = INVEU * downzdnear - MU * _mm512_sqrt_ps(J * J + K * K) * downnear;
					I2 = -MU * (J * X0 + K * Y0);
					Sin2 = _mm512_sincos_ps(&Cos2, I2);
					EQ3 = E2 * Cos2;
					EQ4 = E2 * Sin2;
					res3 += _mm512_reduce_add_ps(EQ3 * W0);
					res4 += _mm512_reduce_add_ps(EQ4 * W0);
				}
				RightTermReal[i] += res1;
				RightTermImag[i] += res2;
				RightTermReal[(2 * NJKBound + 1) * (2 * NJKBound + 1) + i] += res3;
				RightTermImag[(2 * NJKBound + 1) * (2 * NJKBound + 1) + i] += res4;
			}
		}
		else if(tolerance < 0.000001)
		{
			double IntegralTop_X[int(ceil(S * Nw * S * Nw / 8.0)) * 8], IntegralTop_Y[int(ceil(S * Nw * S * Nw / 8.0)) * 8], IntegralTop_W[int(ceil(S * Nw * S * Nw / 8.0)) * 8];
			double IntegralDown_X[int(ceil(S * Nw * S * Nw / 8.0)) * 8], IntegralDown_Y[int(ceil(S * Nw * S * Nw / 8.0)) * 8], IntegralDown_W[int(ceil(S * Nw * S * Nw / 8.0)) * 8];
			double TopNear_F[int(ceil(S * Nw * S * Nw / 8.0)) * 8], TopZDNear_F[int(ceil(S * Nw * S * Nw / 8.0)) * 8], DownNear_F[int(ceil(S * Nw * S * Nw / 8.0)) * 8], DownZDNear_F[int(ceil(S * Nw * S * Nw / 8.0)) * 8];
			for (int i = 0; i < S * Nw * S * Nw; i++)
			{
				int ix = floor(i / (S * Nw));
				int iy = i - ix * (S * Nw);
				IntegralTop_X[i] = IntegralTop[ix][iy][0];
				IntegralTop_Y[i] = IntegralTop[ix][iy][1];
				IntegralTop_W[i] = IntegralTop[ix][iy][3];
				IntegralDown_X[i] = IntegralDown[ix][iy][0];
				IntegralDown_Y[i] = IntegralDown[ix][iy][1];
				IntegralDown_W[i] = IntegralDown[ix][iy][3];
				TopNear_F[i] = TopNear[ix][iy];
				TopZDNear_F[i] = TopZDNear[ix][iy];
				DownNear_F[i] = DownNear[ix][iy];
				DownZDNear_F[i] = DownZDNear[ix][iy];
			}

			for (int i = S * Nw * S * Nw; i < int(ceil(S * Nw * S * Nw / 8.0)) * 8; i++)
			{
				IntegralTop_X[i] = 0.00;
				IntegralTop_Y[i] = 0.00;
				IntegralDown_X[i] = 0.00;
				IntegralDown_Y[i] = 0.00;
				TopNear_F[i] = 0.00;
				TopZDNear_F[i] = 0.00;
				DownNear_F[i] = 0.00;
				DownZDNear_F[i] = 0.00;
				IntegralTop_W[i] = 0.00;
				IntegralDown_W[i] = 0.00;
			}

			double inveu = (1.0 / EU);

			for (int i = min_index; i <= max_index; i++)
			{
				int j = floor(i / ((2 * NJKBound + 1)));
				int k = i - j * ((2 * NJKBound + 1));
				j = j - NJKBound;
				k = k - NJKBound;

				__m512d X0, Y0, W0, topnear, topzdnear, downnear, downzdnear, I1, E1, E2, I2, Sin1, Sin2, Cos1, Cos2, EQ1, EQ2, EQ3, EQ4;
				__m512d MU, J, K, INVEU;

				MU = _mm512_set1_pd(mu);
				J = _mm512_set1_pd(j);
				K = _mm512_set1_pd(k);
				INVEU = _mm512_set1_pd(inveu);

				double res1 = 0.00, res2 = 0.00, res3 = 0.00, res4 = 0.00;

				for (int m = 0; m < S * Nw * S * Nw; m += 8)
				{
					X0 = _mm512_load_pd(&IntegralTop_X[m]);
					Y0 = _mm512_load_pd(&IntegralTop_Y[m]);
					W0 = _mm512_load_pd(&IntegralTop_W[m]);
					topnear = _mm512_load_pd(&TopNear_F[m]);
					topzdnear = _mm512_load_pd(&TopZDNear_F[m]);
					downnear = _mm512_load_pd(&DownNear_F[m]);
					downzdnear = _mm512_load_pd(&DownZDNear_F[m]);

					I1 = -MU * (J * X0 + K * Y0);
					E1 = INVEU * topzdnear + MU * _mm512_sqrt_pd(J * J + K * K) * topnear;
					Sin1 = _mm512_sincos_pd(&Cos1, I1);
					EQ1 = E1 * Cos1;
					EQ2 = E1 * Sin1;
					res1 += _mm512_reduce_add_pd(EQ1 * W0);
					res2 += _mm512_reduce_add_pd(EQ2 * W0);

					E2 = INVEU * downzdnear - MU * _mm512_sqrt_pd(J * J + K * K) * downnear;
					I2 = -MU * (J * X0 + K * Y0);
					Sin2 = _mm512_sincos_pd(&Cos2, I2);
					EQ3 = E2 * Cos2;
					EQ4 = E2 * Sin2;
					res3 += _mm512_reduce_add_pd(EQ3 * W0);
					res4 += _mm512_reduce_add_pd(EQ4 * W0);
				}
				RightTermReal[i] += res1;
				RightTermImag[i] += res2;
				RightTermReal[(2 * NJKBound + 1) * (2 * NJKBound + 1) + i] += res3;
				RightTermImag[(2 * NJKBound + 1) * (2 * NJKBound + 1) + i] += res4;
			}
		}
	}
}

void HSMA2D::SolveLeastSquareProblem(double* C, double** LeftTermReal, double** LeftTermImag, double* RightTermReal, double* RightTermImag, int p, int row)
{

	int Up = floor((p * p - 1) / 2), Down = ceil((p * p - 1) / 2.0);

	int rowAT, columnATA, columnAT;
	double alpha, beta;
	rowAT = Up; columnAT = row; columnATA = Up;

	alpha = 1.0; beta = 0.00;
	double MatrixAT[rowAT * columnAT], MatrixATA[rowAT * columnATA];
	double MatrixAT_New[Down * columnAT], MatrixATA_New[Down * columnATA];

	for (int i = 0; i < rowAT; i++)
		for (int j = 0; j < columnAT; j++)
		{
			MatrixAT[i * columnAT + j] = LeftTermReal[j][i];
		}

	for (int i = 0; i < (rowAT * columnATA); i++) {
		MatrixATA[i] = 0.0;
	}
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowAT, columnATA, columnAT, alpha, MatrixAT, columnAT, MatrixAT, columnAT, beta, MatrixATA, columnATA);
	
	int InfoHelp, VectorHelp[columnATA];

	for (int i = 0; i < columnATA; i++)
		VectorHelp[i] = 0;

	InfoHelp = LAPACKE_dgetrf(CblasRowMajor, columnATA, columnATA, MatrixATA, columnATA, VectorHelp);
	InfoHelp = LAPACKE_dgetri(CblasRowMajor, columnATA, MatrixATA, columnATA, VectorHelp);

	double BB[columnAT], ATB[columnATA], INV_ATA_ATB[columnATA];

	for (int i = 0; i < columnAT; i++)
	{
		BB[i] = RightTermReal[i];
	}
	for (int i = 0; i < columnATA; i++)
	{
		ATB[i] = 0.00;
		INV_ATA_ATB[i] = 0.00;
	}
	cblas_dgemv(CblasRowMajor, CblasNoTrans, rowAT, columnAT, alpha, MatrixAT, columnAT, BB, 1, beta, ATB, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, columnATA, columnATA, alpha, MatrixATA, columnATA, ATB, 1, beta, INV_ATA_ATB, 1);

	C[0] = 0.00;
	for (int i = 0; i < Up; i++)
	{
		C[2 * (i + 1)] = INV_ATA_ATB[i];
	}

	rowAT = Down; columnAT = row; columnATA = Down;

	alpha = 1.0; beta = 0.00;
	

	for (int i = 0; i < rowAT; i++)
		for (int j = 0; j < columnAT; j++)
		{
			MatrixAT_New[i * columnAT + j] = LeftTermImag[j][i];
		}

	for (int i = 0; i < (rowAT * columnATA); i++) {
		MatrixATA_New[i] = 0.0;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowAT, columnATA, columnAT, alpha, MatrixAT_New, columnAT, MatrixAT_New, columnAT, beta, MatrixATA_New, columnATA);

	int VectorHelp_New[columnATA], InfoHelp_New;
	
	for (int i = 0; i < columnATA; i++)
		VectorHelp_New[i] = 0;

	InfoHelp_New = LAPACKE_dgetrf(CblasRowMajor, columnATA, columnATA, MatrixATA_New, columnATA, VectorHelp_New);

	InfoHelp_New = LAPACKE_dgetri(CblasRowMajor, columnATA, MatrixATA_New, columnATA, VectorHelp_New);

	double BB_New[columnAT],ATB_New[columnATA],INV_ATA_ATB_New[columnATA];

	for (int i = 0; i < columnAT; i++)
	{
		BB_New[i] = RightTermImag[i];
	}
	for (int i = 0; i < columnATA; i++)
	{
		ATB_New[i] = 0.00;
		INV_ATA_ATB_New[i] = 0.00;
	}

	cblas_dgemv(CblasRowMajor, CblasNoTrans, rowAT, columnAT, alpha, MatrixAT_New, columnAT, BB_New, 1, beta, ATB_New, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, columnATA, columnATA, alpha, MatrixATA_New, columnATA, ATB_New, 1, beta, INV_ATA_ATB_New, 1);

	for (int i = 0; i < Down; i++)
	{
		C[2 * i + 1] = INV_ATA_ATB_New[i];
	}	
}

double HSMA2D::FinalCalculateEnergyAndForce(double Force[][3], double* Pot, double Source[][3], double* Q, int NSource, double ImageCharge[][5], int ImageNumber, double **Fibonacci, double** QRD, double** QLocalRD, double Gamma, double* C, int p, double Fp, double F, double Rs, double PI, int IF_FMM_FinalPotential, double tolerance)
{
	if (!IF_FMM_FinalPotential)
	{

		if (tolerance > 0.000001)
		{
			float EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
			float EN[NSource], ENX[NSource], ENY[NSource], ENZ[NSource];
			float Q_Image[int(ceil(ImageNumber / 16.0)) * 16];
			for (int j = 0; j < ImageNumber; j++)
			{
				Q_Image[j] = ImageCharge[j][3] * pow(Gamma, abs(ImageCharge[j][4]));
			}
			for (int j = ImageNumber; j<int(ceil(ImageNumber / 16.0)) * 16; j++)
			{
				Q_Image[j] = 0.00;
			}

			double C_New[int(ceil(p * p / 8.0)) * 8];
			for (int i = 0; i<int(ceil(p * p / 8.0)) * 8; i++)
			{
				if (i < p * p)
				{
					C_New[i] = C[i];
				}
				else
				{
					C_New[i] = 0.00;
				}
			}

			float Image_X[int(ceil(ImageNumber / 16.0)) * 16], Image_Y[int(ceil(ImageNumber / 16.0)) * 16], Image_Z[int(ceil(ImageNumber / 16.0)) * 16];
			for (int i = 0; i < int(ceil(ImageNumber / 16.0)) * 16; i++)
			{
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
				}
			}

			#pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = omp_get_num_threads();

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					double QF[int(ceil(p * p/8.0))*8], QFX[int(ceil(p * p / 8.0)) * 8], QFY[int(ceil(p * p / 8.0)) * 8], QFZ[int(ceil(p * p / 8.0)) * 8];
					CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
					EF[i] = 0.00; EFX[i] = 0.00; EFY[i] = 0.00; EFZ[i] = 0.00;
					EN[i] = 0.00; ENX[i] = 0.00; ENY[i] = 0.00; ENZ[i] = 0.00;
					for (int ii = p * p; ii<int(ceil(p * p / 8.0)) * 8; ii++)
					{
						QF[ii] = 0.00;
						QFX[ii] = 0.00;
						QFY[ii] = 0.00;
						QFZ[ii] = 0.00;
					}
					__m512d qf, qfz, qfx, qfy, c;
					for (int j = 0; j < p * p; j = j + 8)
					{
						qf = _mm512_load_pd(&QF[j]);
						qfz = _mm512_load_pd(&QFZ[j]);
						qfx = _mm512_load_pd(&QFX[j]);
						qfy = _mm512_load_pd(&QFY[j]);
						c = _mm512_load_pd(&C_New[j]);
						EF[i] = EF[i] + _mm512_reduce_add_pd(qf * c);
						EFX[i] = EFX[i] + _mm512_reduce_add_pd(qfx * c);
						EFY[i] = EFY[i] + _mm512_reduce_add_pd(qfy * c);
						EFZ[i] = EFZ[i] + _mm512_reduce_add_pd(qfz * c);
					}
					__m512 X0, Y0, Z0, X1, Y1, Z1, Q1, deltax, deltay, deltaz, delta, square, judge, Zero, midterm, delta_square;
					__mmask16 k0;
					X0 = _mm512_set1_ps(Source[i][0]);
					Y0 = _mm512_set1_ps(Source[i][1]);
					Z0 = _mm512_set1_ps(Source[i][2]);
					judge = _mm512_set1_ps(0.000000000001);
					Zero = _mm512_set1_ps(0.00);
					for (int j = 0; j < ImageNumber; j = j + 16)
					{
						X1 = _mm512_load_ps(&Image_X[j]);
						Y1 = _mm512_load_ps(&Image_Y[j]);
						Z1 = _mm512_load_ps(&Image_Z[j]);
						Q1 = _mm512_load_ps(&Q_Image[j]);
						deltax = X1 - X0;
						deltay = Y1 - Y0;
						deltaz = Z1 - Z0;
						square = deltax * deltax + deltay * deltay + deltaz * deltaz;
						k0 = _mm512_cmp_ps_mask(square, judge, _MM_CMPINT_GT);
						delta = _mm512_mask_invsqrt_ps(Zero, k0, square);
						midterm = Q1 * delta;
						delta_square = delta * delta;
						EN[i] += _mm512_reduce_add_ps(midterm);
						ENX[i] += _mm512_reduce_add_ps(midterm * delta_square * deltax);
						ENY[i] += _mm512_reduce_add_ps(midterm * delta_square * deltay);
						ENZ[i] += _mm512_reduce_add_ps(midterm * delta_square * deltaz);

					}
				}

			}

			double Energy = 0.00;
			for (int i = 0; i < NSource; i++)
			{
				Pot[i] = EN[i] + EF[i];
				Energy = Energy + Q[i] * Pot[i];
				Force[i][0] = -(EFX[i] + ENX[i]) * Q[i];
				Force[i][1] = -(EFY[i] + ENY[i]) * Q[i];
				Force[i][2] = -(EFZ[i] + ENZ[i]) * Q[i];
			}

			return Energy;
		}
		else
		{
			double EF[NSource], EFX[NSource], EFY[NSource], EFZ[NSource];
			double EN[NSource], ENX[NSource], ENY[NSource], ENZ[NSource];
			double Q_Image[int(ceil(ImageNumber / 8.0)) * 8];
			
			for (int j = 0; j < ImageNumber; j++)
			{
				Q_Image[j] = ImageCharge[j][3] * pow(Gamma, abs(ImageCharge[j][4]));
			}
			for (int j = ImageNumber; j<int(ceil(ImageNumber /8.0)) * 8; j++)
			{
				Q_Image[j] = 0.00;
			}

			double C_New[int(ceil(p * p / 8.0)) * 8];
			for (int i = 0; i<int(ceil(p * p / 8.0)) * 8; i++)
			{
				if (i < p * p)
				{
					C_New[i] = C[i];
				}
				else
				{
					C_New[i] = 0.00;
				}
			}

			double Image_X[int(ceil(ImageNumber / 8.0)) * 8], Image_Y[int(ceil(ImageNumber / 8.0)) * 8], Image_Z[int(ceil(ImageNumber / 8.0)) * 8];
			for (int i = 0; i < int(ceil(ImageNumber / 8.0)) * 8; i++)
			{	
				if (i < ImageNumber) {
					Image_X[i] = ImageCharge[i][0];
					Image_Y[i] = ImageCharge[i][1];
					Image_Z[i] = ImageCharge[i][2];
				}
				else {
					Image_X[i] = 0.00;
					Image_Y[i] = 0.00;
					Image_Z[i] = 0.00;
				}
			}

		#pragma omp parallel
			{
				int id = omp_get_thread_num();
				int size = omp_get_num_threads();

				int min_atom = id * floor(NSource / size) + 1, max_atom = (id + 1) * floor(NSource / size);
				if (id == size - 1)max_atom = NSource - 1;
				if (id == 0)min_atom = 0;

				for (int i = min_atom; i <= max_atom; i++)
				{
					double QF[int(ceil(p * p / 8.0)) * 8], QFX[int(ceil(p * p / 8.0)) * 8], QFY[int(ceil(p * p / 8.0)) * 8], QFZ[int(ceil(p * p / 8.0)) * 8];
					CalculateMultipoleExpansion(QF, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateZDerivativeMultipoleExpansion(QFZ, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateXDMultipoleExpansion(QFX, p, Source[i][0], Source[i][1], Source[i][2]);
					CalculateYDMultipoleExpansion(QFY, p, Source[i][0], Source[i][1], Source[i][2]);
					EF[i] = 0.00; EFX[i] = 0.00; EFY[i] = 0.00; EFZ[i] = 0.00;
					EN[i] = 0.00; ENX[i] = 0.00; ENY[i] = 0.00; ENZ[i] = 0.00;
					for (int ii = p * p; ii<int(ceil(p * p / 8.0)) * 8; ii++)
					{
						QF[ii] = 0.00;
						QFX[ii] = 0.00;
						QFY[ii] = 0.00;
						QFZ[ii] = 0.00;
					}
					__m512d qf, qfz, qfx, qfy, c;
					for (int j = 0; j < p * p; j = j + 8)
					{
						qf = _mm512_load_pd(&QF[j]);
						qfz = _mm512_load_pd(&QFZ[j]);
						qfx = _mm512_load_pd(&QFX[j]);
						qfy = _mm512_load_pd(&QFY[j]);
						c = _mm512_load_pd(&C_New[j]);
						EF[i] = EF[i] + _mm512_reduce_add_pd(qf * c);
						EFX[i] = EFX[i] + _mm512_reduce_add_pd(qfx * c);
						EFY[i] = EFY[i] + _mm512_reduce_add_pd(qfy * c);
						EFZ[i] = EFZ[i] + _mm512_reduce_add_pd(qfz * c);
					}
					__m512d X0, Y0, Z0, X1, Y1, Z1, Q1, deltax, deltay, deltaz, delta, square, judge, Zero, midterm, delta_square;
					__mmask16 k0;
					X0 = _mm512_set1_pd(Source[i][0]);
					Y0 = _mm512_set1_pd(Source[i][1]);
					Z0 = _mm512_set1_pd(Source[i][2]);
					judge = _mm512_set1_pd(0.000000000001);
					Zero = _mm512_set1_pd(0.00);
					for (int j = 0; j < ImageNumber; j = j + 8)
					{
						X1 = _mm512_load_pd(&Image_X[j]);
						Y1 = _mm512_load_pd(&Image_Y[j]);
						Z1 = _mm512_load_pd(&Image_Z[j]);
						Q1 = _mm512_load_pd(&Q_Image[j]);
						deltax = X1 - X0;
						deltay = Y1 - Y0;
						deltaz = Z1 - Z0;
						square = deltax * deltax + deltay * deltay + deltaz * deltaz;
						k0 = _mm512_cmp_pd_mask(square, judge, _MM_CMPINT_GT);
						delta = _mm512_mask_invsqrt_pd(Zero, k0, square);
						midterm = Q1 * delta;
						delta_square = delta * delta;
						EN[i] += _mm512_reduce_add_pd(midterm);
						ENX[i] += _mm512_reduce_add_pd(midterm * delta_square * deltax);
						ENY[i] += _mm512_reduce_add_pd(midterm * delta_square * deltay);
						ENZ[i] += _mm512_reduce_add_pd(midterm * delta_square * deltaz);
					}
				}
			}

			double Energy = 0.00;
			for (int i = 0; i < NSource; i++)
			{
				Pot[i] = EN[i] + EF[i];
				Energy = Energy + Q[i] * Pot[i];
				Force[i][0] = -(EFX[i] + ENX[i]) * Q[i];
				Force[i][1] = -(EFY[i] + ENY[i]) * Q[i];
				Force[i][2] = -(EFZ[i] + ENZ[i]) * Q[i];
			}

			return Energy;
		}
	}
	else
	{

		double eps = tolerance;
		int ns = ImageNumber;
		int nt = NSource;

		double* source = (double*)malloc(3 * ns * sizeof(double));
		double* target = (double*)malloc(3 * nt * sizeof(double));
		double* charge = (double*)malloc(ns * sizeof(double));

		for (int i = 0; i < ns; i++)
		{
			source[3 * i] = ImageCharge[i][0];
			source[3 * i + 1] = ImageCharge[i][1];
			source[3 * i + 2] = ImageCharge[i][2];
			charge[i] = ImageCharge[i][3] * pow(Gamma, fabs(ImageCharge[i][4]) + 0.00);
		}

		for (int i = 0; i < nt; i++)
		{
			target[3 * i] = Source[i][0];
			target[3 * i + 1] = Source[i][1];
			target[3 * i + 2] = Source[i][2];
		}

		double* pottarg = (double*)malloc(nt * sizeof(double));
		double* gradtarg = (double*)malloc(3 * nt * sizeof(double));

		int ier = 0;
		lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, target, pottarg, gradtarg, &ier);

		int KL = int(2 * F + 2);

		/*             BEGIN HSMA ALGORITHM                          */
		double EF[KL];
		double CenterPara;

		for (int i = 0; i < KL; i++)
		{
			CenterPara = 0.00;
			for (int j = 0; j < p * p; j++)
			{
				CenterPara = CenterPara + (1 / (4 * PI)) * (QRD[i][j] - QLocalRD[i][j]) * C[j];
			}
			EF[i] = CenterPara * Fibonacci[i][3];
		}

		/*            Set FMM parameters          */
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

		lfmm3d_t_c_g_(&eps, &ns, sourceF, chargeF, &nt, target, pottargF, gradtargF, &ier);

		/*	Final Summation	*/
		double Energy = 0.00;;
		for (int i = 0; i < NSource; i++)
		{
			Pot[i] = pottarg[i] + pottargF[i];
			Energy = Energy + Pot[i] * Q[i];
			Force[i][0] = -(gradtarg[3 * i] + gradtargF[3 * i]) * Q[i];
			Force[i][1] = -(gradtarg[3 * i + 1] + gradtargF[3 * i + 1]) * Q[i];
			Force[i][2] = -(gradtarg[3 * i + 2] + gradtargF[3 * i + 2]) * Q[i];
		}

		cout << pottarg[0] << "   " << pottargF[0] << endl;
		return Energy;
	}
}
