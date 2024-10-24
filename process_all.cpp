#include <iostream>
#include <dirent.h>
#include <sys/types.h>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <limits.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TFitResultPtr.h"
#include "TList.h"
#include "TSystem.h"

// #define _draw_waveform
// define _draw_waveform to draw waveforms of both gain channel;

typedef struct
{
    Long64_t CrystalID;
    Double_t Temperature1;
    Double_t Temperature2;
    Double_t LAmplitude[256];
    Double_t HAmplitude[256]; // waveform points
    Double_t LNoise[25];
    Double_t HNoise[25]; // first 25 points of waveform
    Double_t LowGainPedestal;
    Double_t HighGainPedestal; // average of Noise
    Double_t LowGainPeak;
    Double_t HighGainPeak; // peak value of waveform
} CHANNEL_SIGNAL;

typedef struct
{
    Long64_t EventID;
    // Int_t TriggerID; // something goes wrong when reading TriggerID
    Long64_t TimeCode;
    float Time[256];
    float Voltage[5];
    float Current[5];
} CHANNEl_STATUS;

typedef struct
{
    Long64_t CrystalID;
    Double_t Temperature1;
    Double_t Temperature2;
    Double_t LowGainSpectrum;
    Double_t HighGainSpectrum; // peak - pedestal
    Bool_t OverThreshold;      // true if (lowgain_peak > lowgain_pedestal_mean + 10 * lowgain_pedestal_sigma)
} CHANNEL_SPECTRUM;

// ----------------------------------------------------------------------------------
// -------------------- COPY FROM ROOT TUTORIAL '/fit/langaus.C' --------------------
// ------------- Convoluted Landau and Gaussian Fitting Function --------------------
double langaufun(Double_t *x, Double_t *par)
{

    // Fit parameters:
    // par[0]=Width (scale) parameter of Landau density
    // par[1]=Most Probable (MP, location) parameter of Landau density
    // par[2]=Total area (integral -inf to inf, normalization constant)
    // par[3]=Width (sigma) of convoluted Gaussian function
    //
    // In the Landau distribution (represented by the CERNLIB approximation),
    // the maximum is located at x=-0.22278298 with the location parameter=0.
    // This shift is corrected within this function, so that the actual
    // maximum is identical to the MP parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
    Double_t mpshift = -0.22278298;      // Landau maximum location

    // Control constants
    Double_t np = 100.0; // number of convolution steps
    Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow, xupp;
    Double_t step;
    Double_t i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp - xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for (i = 1.0; i <= np / 2; i++)
    {
        xx = xlow + (i - .5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);

        xx = xupp - (i - .5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}

void langaufit(TH1D *his)
{
    // Once again, here are the Landau * Gaussian parameters:
    //   par[0]=Width (scale) parameter of Landau density
    //   par[1]=Most Probable (MP, location) parameter of Landau density
    //   par[2]=Total area (integral -inf to inf, normalization constant)
    //   par[3]=Width (sigma) of convoluted Gaussian function
    //
    // Variables for langaufit call:
    //   his             histogram to fit
    //   startvalues[4]  reasonable start values for the fit
    //   parlimitslo[4]  lower parameter limits
    //   parlimitshi[4]  upper parameter limits

    if (his->GetEntries() < 500)
        return;
    TFitResultPtr fit_res;
    int binMax, binLeft, binRight;
    double thresholdValue, LeftValue, RightValue, LeftleftValue, Rightrightvalue;
    // determine fitting range
    binMax = his->GetMaximumBin();
    thresholdValue = his->GetBinContent(binMax) / 4;
    binLeft = binRight = binMax;
    do
    {
        LeftValue = his->GetBinContent(binLeft);
        RightValue = his->GetBinContent(binRight);
        LeftleftValue = his->GetBinContent(binLeft - 1);
        Rightrightvalue = his->GetBinContent(binRight - 1);
        if (LeftValue > thresholdValue || LeftleftValue > thresholdValue)
            binLeft--;
        if (RightValue > thresholdValue || Rightrightvalue > thresholdValue)
            binRight++;
        if (binLeft == 0 || binRight == his->GetNbinsX())
            break;
    } while (LeftValue > thresholdValue || RightValue > thresholdValue || LeftleftValue > thresholdValue || Rightrightvalue > thresholdValue);
    assert(LeftValue < thresholdValue && RightValue < thresholdValue);
    // determine range done

    // start values, parameter lower limits, parameter upper limits
    double StdDev, MaxPos, Area, Sigma;
    StdDev = his->GetStdDev();
    MaxPos = his->GetMaximumBin() * his->GetBinWidth(0);
    Area = his->Integral("width");
    double startvalues[4] = {StdDev, MaxPos, Area, 1}, parlimitslo[4] = {StdDev * 0.001, MaxPos * 0.1, Area * 0.01, 0.01}, parlimitshi[4] = {StdDev, MaxPos * 10, Area * 100, 100};

    TString FunName = "lan_gaus_conv";
    TF1 *ffitold = (TF1 *)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold)
        delete ffitold;

    TF1 *ffit = new TF1(FunName.Data(), langaufun, his->GetBinCenter(binLeft), his->GetBinCenter(binRight), 4);
    ffit->SetParameters(startvalues);
    ffit->SetParNames("Width", "MP", "Area", "GSigma");

    for (int i = 0; i < 4; i++)
    {
        ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
    }

    fit_res = his->Fit(FunName, "QS RB"); // fit within specified range, use ParLimits
}
// ----------------------------------------------------------------------------------

// input a histogram, search for its 1/4 maximum width as fitting window.
void landau_fit(TH1D *);
// input a graph of high/low ratio, fit with linear function from 0 to 0.9*maximum.
void linear_fit(TGraph *);
// process Step1.root to generate Step1_processed.root
/*
CH1-5 CH2-5 ***** ***** CH5-5             CH05  CH10  ****  ****  CH25
CH1-4 ***** ***** ***** *****             CH04  ****  ****  ****  ****
CH1-3 ***** ***** ***** *****   ====>>    CH03  ****  ****  ****  ****
CH1-2 CH2-2 ***** ***** *****             CH02  CH07  ****  ****  ****
CH1-1 CH2-1 CH3-1 ***** CH5-1             CH01  CH06  CH11  CH16  CH21
*/
void process_data(TString);
// input Step1_processed.root, fit 25 channels with Landau function;
void peak_fit(TString);
// draw spectra of 25 channels in *_processeed.root file
void draw_spectra(TString);
// High/Low gain ratio for 25 channels (subtract pedestal)
// input *_processed.root file
void gain_ratio(TString);

// process all step1 files in a directory to generated corresponding *_processed.root files, and hadd to 'merged.root' file.
// parameter is the path of directory either absolute or relative one, default is current work directory.
void process_all(TString path_name = "./")
{
    gStyle->SetOptFit(1111);
    TFile *infile;
    TString filename;
    TString file_path;
    DIR *dirp = NULL;
    struct dirent *dir_entry = NULL;
    char current_path[PATH_MAX] = "";
    if ((dirp = opendir(path_name.Data())) == NULL)
    {
        std::cout << "convert to absolute path" << std::endl;
        getcwd(current_path, PATH_MAX);
        path_name = Form("%s/%s", current_path, path_name.Data());
        if ((dirp = opendir(path_name.Data())) == NULL)
        {
            std::cerr << "find path error!" << std::endl;
            return;
        }
    }
    while ((dir_entry = readdir(dirp)) != NULL)
    {
        filename = dir_entry->d_name;
        if (filename.Contains("_processed") || (!filename.Contains(".root")))
            continue;
        else if (filename == "merged.root")
            gSystem->Exec(Form("rm %s/merged.root", path_name.Data()));
        else
        {
            std::cout << "reading file: " << filename << std::endl;
            file_path = Form("%s/%s", path_name.Data(), filename.Data());
            process_data(file_path);
            peak_fit(Form("%s_processed.root", file_path.Remove(file_path.Last('.')).Data()));
        }
    }
    gSystem->Exec(Form("hadd %s/merged.root %s/*_processed.root", path_name.Data(), path_name.Data()));
    peak_fit(Form("%s/merged.root", path_name.Data()));
}

void landau_fit(TH1D *his)
{
    if (his->GetEntries() < 500)
        return;
    TFitResultPtr fit_res;
    int binMax, binLeft, binRight;
    double thresholdValue, LeftValue, RightValue, LeftleftValue, Rightrightvalue;
    binMax = his->GetMaximumBin();
    thresholdValue = his->GetBinContent(binMax) / 4;
    binLeft = binRight = binMax;
    do
    {
        LeftValue = his->GetBinContent(binLeft);
        RightValue = his->GetBinContent(binRight);
        LeftleftValue = his->GetBinContent(binLeft - 1);
        Rightrightvalue = his->GetBinContent(binRight - 1);
        if (LeftValue > thresholdValue || LeftleftValue > thresholdValue)
            binLeft--;
        if (RightValue > thresholdValue || Rightrightvalue > thresholdValue)
            binRight++;
        if (binLeft == 0 || binRight == his->GetNbinsX())
            break;
    } while (LeftValue > thresholdValue || RightValue > thresholdValue || LeftleftValue > thresholdValue || Rightrightvalue > thresholdValue);
    assert(LeftValue < thresholdValue && RightValue < thresholdValue);
    fit_res = his->Fit("landau", "QS", "", his->GetBinCenter(binLeft), his->GetBinCenter(binRight));
}

void linear_fit(TGraph *gr)
{
    if (gr->GetN() < 500)
        return;
    gr->Sort();
    TFitResultPtr fit_res;
    double Xleft = 0, threshold;
    int eventNum, thresholdID = 0;
    eventNum = gr->GetN();
    threshold = gr->GetPointY(gr->GetN() - 1) * 0.95;
    for (int i = 0; i < gr->GetN(); i++)
    {
        if (gr->GetPointY(i) > threshold)
        {
            if (gr->GetPointY(i + 1) > threshold)
            {
                thresholdID = i;
                Xleft = gr->GetPointX(i);
                break;
            }
        }
    }
    assert(!(thresholdID == 0 || Xleft == 0));
    fit_res = gr->Fit("pol1", "QS", "", 0, Xleft);
}

void process_data(TString filename)
{
    gStyle->SetOptFit(1111);
    TFile *infile = TFile::Open(filename.Data(), "READ");
    TString outfile_name = Form("%s_processed.root", filename.Remove(filename.Last('.')).Data());
    TFile *outfile = TFile::Open(outfile_name.Data(), "RECREATE");
    if (!infile->IsOpen() && !outfile->IsOpen())
        std::cerr << "Open file error!" << std::endl;
    std::cout << "generating processed file: " << outfile_name << std::endl;
    TTree *tr = (TTree *)infile->Get("decode_data");
    CHANNEL_SIGNAL channel_signal[25];
    CHANNEl_STATUS channel_status;
    CHANNEL_SPECTRUM channel_spectrum[25];
    int EventNum = tr->GetEntries();
    tr->SetBranchAddress("EventID", &channel_status.EventID);
    // tr->SetBranchAddress("TriggerID", &channel_status.TriggerID);
    tr->SetBranchAddress("TimeCode", &channel_status.TimeCode);
    tr->SetBranchAddress("Time", channel_status.Time);
    for (int i = 0; i < 5; i++)
    {
        tr->SetBranchAddress(Form("Voltage_%i", i + 1), &channel_status.Voltage[i]);
        tr->SetBranchAddress(Form("Current_%i", i + 1), &channel_status.Current[i]);
    }
    TTree *tr_pro = new TTree("processed_data", "processed_data");
    tr_pro->SetDirectory(outfile);
    TString leaflist_status = "EventID/L:TimeCode/L:Time[256]/F:Voltage[5]/F:Current[5]/F";
    tr_pro->Branch("Channel_Status", &channel_status, leaflist_status.Data());
    TString channel_name;
    TString leaflist_signal = "";
    leaflist_signal += "CrystalID/L";
    leaflist_signal += ":Temperature1/D";
    leaflist_signal += ":Temperature2/D";
    leaflist_signal += ":LowGainSpectrum/D";
    leaflist_signal += ":HighGainSpectrum/D";
    leaflist_signal += ":OverThreshold/B";
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
            channel_name = Form("Hit_%i_%i", i + 1, j + 1);
            tr->SetBranchAddress(channel_name.Data(), channel_signal + 5 * i + j);
            tr_pro->Branch(Form("Channel_Spectrum_%i", 5 * i + j + 1), channel_spectrum + 5 * i + j, leaflist_signal);
        }

#ifdef _draw_waveform
    TCanvas *can = new TCanvas("can", "can", 900, 600);
    TGraph *gr_hg = new TGraph(256);
    TGraph *gr_lg = new TGraph(256);
    gr_hg->SetTitle("HG waveform;time/ns;ADC value");
    gr_lg->SetTitle("LG waveform;time/ns;ADC value");
    double time[256];
    for (int i = 0; i < 256; i++)
        time[i] = i * 8;
#endif

    TFitResultPtr fit_res;
    TH1D *LowGain_pedestal[25], *HighGain_pedestal[25];
    for (int i = 0; i < 25; i++)
    {
        LowGain_pedestal[i] = new TH1D();
        LowGain_pedestal[i]->SetDirectory(NULL);
        LowGain_pedestal[i]->SetNameTitle(Form("LowGain_pedestal_%i", i + 1), Form("Channel_%i", i + 1));
        LowGain_pedestal[i]->SetBins(300, 2200, 2500);
        HighGain_pedestal[i] = new TH1D();
        HighGain_pedestal[i]->SetDirectory(NULL);
        HighGain_pedestal[i]->SetNameTitle(Form("HighGain_pedestal_%i", i + 1), Form("Channel_%i", i + 1));
        HighGain_pedestal[i]->SetBins(600, 2000, 2600);
    }

    for (int i = 0; i < EventNum; i++)
    {
        if (i % (EventNum / 10) == 0 && i != 0)
            std::cout << "Process: " << 5 * i / (EventNum / 10) << "%" << std::endl;
        tr->GetEntry(i);

#ifdef _draw_waveform
        // draw channel 13
        if (channel_signal[13].LowGainPeak > 3000)
        {
            for (int j = 0; j < 256; j++)
            {
                gr_hg->SetPoint(j, time[j], channel_signal[13].HAmplitude[j]);
            }
            gr_hg->Draw();
            can->Update();
            can->SaveAs(Form("highgain_waveform_%i.png", i));
            can->Clear();
            sleep(2);
            for (int j = 0; j < 256; j++)
            {
                gr_lg->SetPoint(j, time[j], channel_signal[13].LAmplitude[j]);
            }
            gr_lg->Draw();
            can->Update();
            can->SaveAs(Form("lowgain_waveform_%i.png", i));
            sleep(2);
            can->Clear();
        }
#endif
        // fill pedestal of each channels
        for (int j = 0; j < 25; j++)
        {
            for (int k = 0; k < 25; k++)
            {
                LowGain_pedestal[j]->Fill(channel_signal[j].LAmplitude[k]);
                HighGain_pedestal[j]->Fill(channel_signal[j].HAmplitude[k]);
            }
        }
    }
    // fit to get pedestal mean and sigma, set threshold = 10 * sigma;
    double LowGain_mean[25], LowGain_sigma[25], LowGain_threshold[25];
    double HighGain_mean[25], HighGain_sigma[25], HighGain_threshold[25];
    for (int i = 0; i < 25; i++)
    {
        fit_res = LowGain_pedestal[i]->Fit("gaus", "QS0");
        LowGain_mean[i] = fit_res->Parameter(1);
        LowGain_sigma[i] = fit_res->Parameter(2);
        LowGain_threshold[i] = LowGain_sigma[i] * 10;
        fit_res = HighGain_pedestal[i]->Fit("gaus", "QS0");
        HighGain_mean[i] = fit_res->Parameter(1);
        HighGain_sigma[i] = fit_res->Parameter(2);
        HighGain_threshold[i] = HighGain_sigma[i] * 10;
    }
    // TCanvas *can = new TCanvas("can", "can", 1600, 900);
    // can->Divide(5, 5);
    // for (int i = 0; i < 5; i++)
    // {
    //     for (int j = 0; j < 5; j++)
    //     {
    //         can->cd(5 * i + j + 1);
    //         LowGain_pedestal[5 * j + 4 - i]->Draw();
    //         can->Update();
    //         std::cout << LowGain_mean[5 * j + 4 - i] << '\t';
    //     }
    //     std::cout << std::endl;
    // }

    // calculate spectrum = peak - pedestal_mean, and fill tr_pro
    for (int i = 0; i < EventNum; i++)
    {
        if (i % (EventNum / 10) == 0 && i != 0)
            std::cout << "Process: " << 5 * (i + EventNum) / (EventNum / 10) << "%" << std::endl;
        tr->GetEntry(i);
        for (int j = 0; j < 25; j++)
        {
            channel_spectrum[j].CrystalID = channel_signal[j].CrystalID;
            channel_spectrum[j].Temperature1 = channel_signal[j].Temperature1;
            channel_spectrum[j].Temperature2 = channel_signal[j].Temperature2;
            channel_spectrum[j].HighGainSpectrum = channel_signal[j].HighGainPeak - HighGain_mean[j];
            channel_spectrum[j].LowGainSpectrum = channel_signal[j].LowGainPeak - LowGain_mean[j];
            if (channel_spectrum[j].LowGainSpectrum > LowGain_threshold[j])
                channel_spectrum[j].OverThreshold = true;
            else
                channel_spectrum[j].OverThreshold = false;
        }
        tr_pro->Fill();
    }
    infile->Close();
    outfile->Write();
    outfile->Close();
}

void peak_fit(TString filename)
{
    gStyle->SetOptFit(1111);
    TFile *infile = TFile::Open(filename.Data(), "UPDATE");
    if (!infile->IsOpen())
        std::cerr << "open file error!" << std::endl;
    std::cout << "fitting histograms:" << std::endl;
    TTree *tr = (TTree *)infile->Get("processed_data");
    CHANNEL_SPECTRUM channel_spectrum[25];
    for (int i = 0; i < 25; i++)
    {
        tr->SetBranchAddress(Form("Channel_Spectrum_%i", i + 1), &channel_spectrum[i]);
    }
    TH1D *LowEnergy_his[25];
    for (int i = 0; i < 25; i++)
    {
        LowEnergy_his[i] = new TH1D();
        LowEnergy_his[i]->SetDirectory(infile);
        LowEnergy_his[i]->SetNameTitle(Form("lowspectrum_%i", i + 1), Form("channel_%i_spectrum", i + 1));
        LowEnergy_his[i]->SetBins(500, 0, 2000);
    }
    int EventNum = tr->GetEntries();
    int judge;
    for (int i = 0; i < EventNum; i++)
    {
        if (i % (EventNum / 20) == 0 && i != 0)
            std::cout << "Process: " << 100 * i / EventNum << "%" << std::endl;
        tr->GetEntry(i);
        judge = 0;
        for (int j = 0; j < 25; j++)
            judge += channel_spectrum[j].OverThreshold;
        if (judge == 1) // select single hit events
        {
            for (int j = 0; j < 25; j++)
            {
                if (channel_spectrum[j].OverThreshold) // draw only hitted channel
                    LowEnergy_his[j]->Fill(channel_spectrum[j].LowGainSpectrum);
            }
        }
    }
    double mean[25] = {0};
    TF1 *flandau;
    for (int i = 0; i < 25; i++)
    {
        // landau_fit(LowEnergy_his[i]);
        langaufit(LowEnergy_his[i]);
        if ((flandau = (TF1 *)((TList *)LowEnergy_his[i]->GetListOfFunctions())->FindObject("landau")) || (flandau = (TF1 *)((TList *)LowEnergy_his[i]->GetListOfFunctions())->FindObject("lan_gaus_conv")))
            mean[i] = flandau->GetParameter(1);
    }
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            std::cout << mean[5 * j + 4 - i] << '\t';
        }
        std::cout << std::endl;
    }
    infile->Write(0, TObject::kOverwrite);
    infile->Close();
}

void draw_spectra(TString filename)
{
    TFile *infile = TFile::Open(filename.Data(), "read");
    TH1D *his[25];
    double mean[25] = {0};
    TF1 *flandau;
    for (int i = 0; i < 25; i++)
    {
        his[i] = (TH1D *)infile->Get(Form("lowspectrum_%i", i + 1));
        his[i]->GetXaxis()->SetTitle("QDC_channel");
        his[i]->GetYaxis()->SetTitle("count");
        if ((flandau = (TF1 *)((TList *)his[i]->GetListOfFunctions())->FindObject("landau")) || (flandau = (TF1 *)((TList *)his[i]->GetListOfFunctions())->FindObject("lan_gaus_conv")))
            mean[i] = flandau->GetParameter(1);
    }
    TCanvas *can = new TCanvas("can", "can", 1600, 900);
    can->Divide(5, 5);
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            can->cd(5 * i + j + 1);
            his[5 * j + 4 - i]->Draw();
            can->Update();
            std::cout << mean[5 * j + 4 - i] << '\t';
        }
        std::cout << std::endl;
    }
    // infile->Close();
}

void gain_ratio(TString filename)
{
    gStyle->SetOptFit(1111);
    TFile *infile = TFile::Open(filename.Data(), "READ");
    if (!infile->IsOpen())
        std::cerr << "open file error!" << std::endl;
    TTree *tr = (TTree *)infile->Get("processed_data");
    CHANNEL_SPECTRUM channel_spectrum[25];
    TGraph *ratio_gr[25];
    TH2D *ratio_his[25];
    int EventNum = tr->GetEntries();
    for (int i = 0; i < 25; i++)
    {
        tr->SetBranchAddress(Form("Channel_Spectrum_%i", i + 1), &channel_spectrum[i]);
        ratio_gr[i] = new TGraph(EventNum);
        ratio_gr[i]->SetTitle(Form("channel_%i;Highgain;Lowgain", i + 1));
    }
    for (int i = 0; i < EventNum; i++)
    {
        if (i % (EventNum / 20) == 0 && i != 0)
            std::cout << "Process: " << 100 * i / EventNum << "%" << std::endl;
        tr->GetEntry(i);
        for (int j = 0; j < 25; j++)
        {
            ratio_gr[j]->SetPoint(i, channel_spectrum[j].LowGainSpectrum, channel_spectrum[j].HighGainSpectrum);
        }
    }
    double slope[25], intercept[25];
    for (int i = 0; i < 25; i++)
    {
        linear_fit(ratio_gr[i]);
        if (TF1 *func = (TF1 *)ratio_gr[i]->FindObject("pol1"))
        {
            intercept[i] = func->GetParameter(0);
            slope[i] = func->GetParameter(1);
        }
    }
    TCanvas *can = new TCanvas("can", "can", 1600, 900);
    can->Divide(5, 5);
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            can->cd(5 * i + j + 1);
            ratio_gr[5 * j + 4 - i]->Draw("ap");
            can->Update();
            std::cout << slope[5 * j + 4 - i] << '\t';
        }
        std::cout << std::endl;
    }
    infile->Close();
}
