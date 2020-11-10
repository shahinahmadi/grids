#include <iostream>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <TGraphErrors.h>
#include <TLatex.h>

double time_average(double* vals) {
	return (vals[0] + vals[1]+vals[2]) / 2;
}

double time_error(double* vals) { 
	double x1 = TMath::Abs(vals[0] - vals[1])/2;

	return x1;
}

double ti2secs(double min, double sec) {
	return min*60 + sec;
}

double rateTotalFitFnc(double* x, double* p) {
	return p[0]*TMath::Exp(-p[1]*x[0]);
}

void Analysis() {

	   	TCanvas *c1 = new TCanvas("c1","the fit canvas", 800, 600);
   c1->SetFillColor(33);
   c1->SetFrameFillColor(255);
   c1->SetGrid();

	const int n = 6;
	double counts[n] = { 193, 155, 141, 123, 105, 106 };
	// double counts[n] = {1341, 776, 3793.0, 105, 96, 35};
	double delay[n] = {99.7, 150.5, 200.0, 250.3, 299.7, 349.7};
	double delay_err[n] = { 0.098, 0.117, 0.121, 0.127, 0.103, 0.129 };
	// double delay[n] = {50.177, 100.3, 150.4, 200.0, 250.9, 299.7};
	double window_time[n] = { 24.81, 25.62, 25.24, 25.20, 25.24, 25.58 };
	double window_time_err[n] = { 0.169, 0.170, 0.090, 0.180, 0.170, 0.180 };
	// double delay_err[n] = { 24.5, 24.50, 24.5, 24.7, 24.5, 25.65 };


	double time_data[n][3] =  { 
	// 	{ti2secs(5, 1.28), ti2secs(5, 1.46), ti2secs(5, 1.86)},
	// 	{ti2secs(5, 8.13), ti2secs(5, 7.47), ti2secs(5, 7.2)},
	// 	{ti2secs(60+8, 57.65), ti2secs(60+8, 57.46), ti2secs(60+8, 58.03)}, 
	// 	{ti2secs(5, 00.88), ti2secs(5, 0.87), ti2secs(5, 0.27)},
	// 	{ti2secs(5, 33.95), ti2secs(5, 33.95), ti2secs(5, 33.06)}, 
	// 	{ti2secs(4, 59.52), ti2secs(5, 00.93), ti2secs(3, 0.9)}
	// };
	{ti2secs(5, 28.13), ti2secs(5, 27.90), ti2secs(5, 27.90)},
	{ti2secs(5, 09.96), ti2secs(5, 09.55), ti2secs(5, 27.90)}, 
	{ti2secs(5, 14.00), ti2secs(5, 24.20), ti2secs(5, 27.90)},
	{ti2secs(5, 20.05), ti2secs(5, 20.47), ti2secs(5, 27.90)}, 
	{ti2secs(5, 13.64), ti2secs(5, 23.42), ti2secs(5, 27.90)},
	{ti2secs(5, 08.21), ti2secs(5, 08.60), ti2secs(5, 27.90)}};

	double* time_avg = new double[n];
	double* time_err = new double[n];
	int i = 0;
	for(auto time : time_data) {
		time_avg[i] = time_average(time);
		time_err[i] = time_error(time);
		i++;
	}

	double* count_err = new double[n];
	i = 0;
	for(auto count: counts) {
		double r_err_t = time_err[i] / time_avg[i];
		double std = (1/count) + r_err_t*r_err_t;
		counts[i] = count / time_avg[i];
		count_err[i] = TMath::Sqrt(std)*counts[i]; 

		std::cout << time_err[i] << std::endl;

		i++;
	}

	TGraphErrors* plot = new TGraphErrors(n, delay, counts, nullptr, count_err);
	plot->SetTitle("BiPo Decay;Time delay (ns);Counts per second");
	plot->SetMarkerStyle(21);
	plot->SetMarkerSize(0.8);


   	auto FitFcn =  new TF1("fitFunction", rateTotalFitFnc, 0, 350, 2);

   	// FitFcn->FixParameter(0, 0.1);
   	// FitFcn->FixParameter(1, 0.1);


	FitFcn->SetParameters(0.1, 0.1);
	FitFcn->SetParNames("A", "L Noise");
	plot->Fit("fitFunction", "M");

   	plot->DrawClone("APE");
	FitFcn->DrawClone("Same");

   	c1->Draw();


	delete [] time_avg, time_err;


   	return c1;
} 
