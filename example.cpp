#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <chrono>

//#include "TFile.h"
#include "TString.h"

#include "gsl/gsl_math.h"

#include "BH_Deu.h"

using namespace std;

int main()
{
	auto starttime = chrono::system_clock::now();
	auto endtime = starttime;
	auto durtime = chrono::duration_cast<chrono::minutes>(endtime - starttime);

	TString filename;
	int Nfiles = 500;
	int j0 = 0;

	// electron: lepType -> 0      muon: lepType -> 1
	int lepType = 0;
	Long64_t simNum = 10000000;

	int physEveNum = 0;
	double var[10];
	double weight = 0.0;
	TLorentzVector pgamFv, pelecFV;
	TLorentzVector pout[3];
	TLorentzVector kf[2];

	double Ebeam = 8.5;

	// initialize the model for Bethe-Heitler process
	BH_deuteron::setModel(lepType, Ebeam);

	// for Bremsstrahlung Photon
	double kmin = 7.2;
	double kmax = Ebeam;
	
	incidentPhoton::SetBremsstrahlung();

	FILE* f;
	TString name;

	filename = "BremPhoton/Event_BremPhoton_";

	for (int j = j0; j < j0 + Nfiles; j++) 
	{
		cout << j << " / " << "[" << j0 << "," << j0 + Nfiles << "]" << endl;

		name = filename + Form("%.4d.dat", j);
		f = fopen(name.Data(), "w");
		fprintf(f, "Ebeam=%.2fGeV,Nsim= %lld\n", Ebeam, simNum);

		for (Long64_t i = 0; i < simNum; i++) 
		{
			// 15cm LD2 target,   0.01 factor is included in the distribution function
			weight = incidentPhoton::BremsstrahlungPhoton(&pgamFv, kmin, kmax, Ebeam) * 1.95 / 2.0;
			weight *= BH_deuteron::GetBHdeu(&pgamFv, pout, var);

			if (weight > 0.0) 
			{
				physEveNum += 1;

				fprintf(f, "%.6E\n", weight);
				fprintf(f, "q:\t%.6E\t%.6E\t%.6E\t%.6E\n", pgamFv.X(), pgamFv.Y(), pgamFv.Z(), pgamFv.E());
				fprintf(f, "e+:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[2].X(), pout[2].Y(), pout[2].Z(), pout[2].E());
				fprintf(f, "e-:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[1].X(), pout[1].Y(), pout[1].Z(), pout[1].E());
				fprintf(f, "p:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[0].X(), pout[0].Y(), pout[0].Z(), pout[0].E());
			}
		}
		fclose(f);
		cout << physEveNum << endl;
	}
	cout << "the physical point number for Bremsstrahlung photon is " << physEveNum << endl;


	// for electroproduction
	double degtorad = M_PI / 180.0;

	filename = "ElecPhoton/Event_Electro_";
	incidentPhoton::perange[0] = 0.0;
	incidentPhoton::perange[1] = Ebeam - 7.2;

	pelecFV.SetXYZM(0., 0., Ebeam, leptonTensor::melec);

	physEveNum = 0;

	for (int j = j0; j < j0 + Nfiles; j++)
	{
		cout << j << " / " << "[" << j0 << "," << j0 + Nfiles << "]" << endl;

		name = filename + Form("%.4d.dat", j);
		f = fopen(name.Data(), "w");
		fprintf(f, "Ebeam=%.2fGeV,Nsim= %lld\n", Ebeam, simNum);

		incidentPhoton::cthrange[0] = cos(5.0 * degtorad);
		incidentPhoton::cthrange[1] = cos(0.0 * degtorad);
		for (Long64_t i = 0; i < simNum; i++)
		{
			// 15cm LD2 target,   0.01 factor is included in the distribution function
			weight = incidentPhoton::VirtualPhoton(&pelecFV, kf);
			pgamFv.SetPxPyPzE(kf[1].Px(), kf[1].Py(), kf[1].Pz(), kf[1].E());
			weight *= BH_deuteron::GetBHdeu(&pgamFv, pout, var);

			if (weight > 0.0)
			{
				physEveNum += 1;

				fprintf(f, "%.6E\n", weight);
				fprintf(f, "e':\t%.6E\t%.6E\t%.6E\t%.6E\n", kf[0].X(), kf[0].Y(), kf[0].Z(), kf[0].E());
				fprintf(f, "e+:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[2].X(), pout[2].Y(), pout[2].Z(), pout[2].E());
				fprintf(f, "e-:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[1].X(), pout[1].Y(), pout[1].Z(), pout[1].E());
				fprintf(f, "p:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[0].X(), pout[0].Y(), pout[0].Z(), pout[0].E());
			}
		}

		incidentPhoton::cthrange[0] = cos(10.0 * degtorad);
		incidentPhoton::cthrange[1] = cos(5.0 * degtorad);
		for (Long64_t i = 0; i < simNum; i++)
		{
			weight = incidentPhoton::VirtualPhoton(&pelecFV, kf);
			pgamFv.SetPxPyPzE(kf[1].Px(), kf[1].Py(), kf[1].Pz(), kf[1].E());
			weight *= BH_deuteron::GetBHdeu(&pgamFv, pout, var);

			if (weight > 0.0)
			{
				physEveNum += 1;

				fprintf(f, "%.6E\n", weight);
				fprintf(f, "e':\t%.6E\t%.6E\t%.6E\t%.6E\n", kf[0].X(), kf[0].Y(), kf[0].Z(), kf[0].E());
				fprintf(f, "e+:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[2].X(), pout[2].Y(), pout[2].Z(), pout[2].E());
				fprintf(f, "e-:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[1].X(), pout[1].Y(), pout[1].Z(), pout[1].E());
				fprintf(f, "p:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[0].X(), pout[0].Y(), pout[0].Z(), pout[0].E());
			}
		}

		incidentPhoton::cthrange[0] = cos(30.0 * degtorad);
		incidentPhoton::cthrange[1] = cos(10.0 * degtorad);
		for (Long64_t i = 0; i < simNum; i++)
		{
			weight = incidentPhoton::VirtualPhoton(&pelecFV, kf);
			pgamFv.SetPxPyPzE(kf[1].Px(), kf[1].Py(), kf[1].Pz(), kf[1].E());
			weight *= BH_deuteron::GetBHdeu(&pgamFv, pout, var);

			if (weight > 0.0)
			{
				physEveNum += 1;

				fprintf(f, "%.6E\n", weight);
				fprintf(f, "e':\t%.6E\t%.6E\t%.6E\t%.6E\n", kf[0].X(), kf[0].Y(), kf[0].Z(), kf[0].E());
				fprintf(f, "e+:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[2].X(), pout[2].Y(), pout[2].Z(), pout[2].E());
				fprintf(f, "e-:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[1].X(), pout[1].Y(), pout[1].Z(), pout[1].E());
				fprintf(f, "p:\t%.6E\t%.6E\t%.6E\t%.6E\n", pout[0].X(), pout[0].Y(), pout[0].Z(), pout[0].E());
			}
		}
		fclose(f);
	}
	cout << "the physical point number for electroproduction photon is " << physEveNum << endl;

	endtime = chrono::system_clock::now();
	durtime = chrono::duration_cast<chrono::minutes>(endtime - starttime);
	cout << "the execution time is " << durtime.count() << " minutes " << endl;

	return 0;
}