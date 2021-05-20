#define lqtoxe_cxx
#include "lqtoxe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

const int arraysize = 8;
float xmax[arraysize];
float xmin[arraysize];
void lqtoxe::Loop()
{
    int massPoint;
    cout << "massPoint=" << endl;       //enter masspoint
    cin >> massPoint;
    cout << "massPoint=" << massPoint << endl;
    
    TH1F *histo[arraysize];
    
    string hist[arraysize] = {"hist_Ele1_Pt","hist_Ele1_Pt_Ele35","hist_Ele1_Pt_Ele35/Photon200","hist_Ele1_Pt_Ele35/Photon200/Ele115","hist_Ele2_Pt","hist_Ele2_Pt_Ele35","hist_Ele2_Pt_Ele35/Photon200","hist_Ele2_Pt_Ele35/Photon200/Ele115"};
    
    string xtitiles[arraysize] = {"Ele1_Pt","Ele1_Pt_Ele35","Ele1_Pt_Ele35/Photon200","Ele1_Pt_Ele35/Photon200/Ele115","Ele2_Pt","Ele2_Pt_Ele35","Ele2_Pt_Ele35/Photon200","Ele2_Pt_Ele35/Photon200/Ele115"};
    
    string titiles[arraysize]={"Ele1_Pt","Ele1_Pt_Ele35","Ele1_Pt_Ele35/Photon200","Ele1_Pt_Ele35/Photon200/Ele115","Ele2_Pt","Ele2_Pt_Ele35","Ele2_Pt_Ele35/Photon200","Ele2_Pt_Ele35/Photon200/Ele115"};
    
    if(massPoint==500) {
        xmin[0]=0; xmin[1]=0; xmin[2]=0; xmin[3]=0; xmin[4]=1; xmin[5]=1; xmin[6]=1; xmin[7]=1;
    }
    else if(massPoint==1000) {
        xmin[0]=0; xmin[1]=0; xmin[2]=0; xmin[3]=0; xmin[4]=1; xmin[5]=1; xmin[6]=1; xmin[7]=1;
    }
    else if(massPoint==1500) {
        xmin[0]=0; xmin[1]=0; xmin[2]=0; xmin[3]=0; xmin[4]=1; xmin[5]=1; xmin[6]=1; xmin[7]=1;
    }
    else if(massPoint==2000) {
        xmin[0]=0; xmin[1]=0; xmin[2]=0; xmin[3]=0; xmin[4]=1; xmin[5]=1; xmin[6]=1; xmin[7]=1;
    }
    else {
        cout << "wrong mass given" << endl;
    }
    
    
    if(massPoint==500) {
        xmax[0]=2000; xmax[1]=2000; xmax[2]=2000; xmax[3]=2000; xmax[4]=2000; xmax[5]=2000; xmax[6]=2000; xmax[7]=2000;
    }
    else if(massPoint==1000) {
        xmax[0]=2000; xmax[1]=2000; xmax[2]=2000; xmax[3]=2000; xmax[4]=2000; xmax[5]=2000; xmax[6]=2000; xmax[7]=2000;
    }
    else if(massPoint==1500) {
        xmax[0]=2000; xmax[1]=2000; xmax[2]=2000; xmax[3]=2000; xmax[4]=2000; xmax[5]=2000; xmax[6]=2000; xmax[7]=2000;
    }
    else if(massPoint==2000) {
        xmax[0]=3000; xmax[1]=3000; xmax[2]=3000; xmax[3]=3000; xmax[4]=3000; xmax[5]=3000; xmax[6]=3000; xmax[7]=3000;
    }
    else {
        cout << "wrong mass given" << endl;
    }
    
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
    
    for(int i = 0; i<arraysize; i++)
    {
        histo[i] = new TH1F(hist[i].c_str(),xtitiles[i].c_str(),100,xmin[i],xmax[i]);
        histo[i]->Sumw2();
    }
    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry > 47000) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       float weight=Weight*puWeight;
       
       if (Ele2_Pt > 35){
           histo[0]->Fill(Ele1_Pt,weight);
           histo[4]->Fill(Ele2_Pt,weight);
       
       if(H_Ele35_WPTight==1){
           histo[1]->Fill(Ele1_Pt,weight);
           histo[5]->Fill(Ele2_Pt,weight);
       }
       if(H_Ele35_WPTight==1 || H_Photon200==1){
           histo[2]->Fill(Ele1_Pt,weight);
           histo[6]->Fill(Ele2_Pt,weight);
       }
       if(H_Ele35_WPTight==1 || H_Photon200==1 || H_Ele115_CIdVT_GsfIdT==1){
           histo[3]->Fill(Ele1_Pt,weight);
           histo[7]->Fill(Ele2_Pt,weight);
       }
   }
   }

    
    
    TFile lqtoxe("lqtoxe.root","recreate");
    lqtoxe.cd();

   
    for(int i = 0; i<arraysize; i++){
        histo[i]->Write();
    }
    
    
    
    
    
    TCanvas* can = new TCanvas(("can"+hist[0]).c_str(),("can"+hist[0]).c_str());
    can->Divide(2,1);   //divide canvas to two part(left and right)
    can->cd(1);   //go to canvas1
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.1); // Upper and lower plot are not joined
    pad1->SetGridx();         // Vertical grid
  //  pad1->SetLogy();
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
   
    
    histo[0]->SetLineColor(1);
    histo[0]->SetLineWidth(2);
    histo[0]->SetName("Ele1_Pt");         //set name
    histo[0]->GetXaxis()->SetTitle(xtitiles[0].c_str());    //set xaxis
    histo[0]->GetYaxis()->SetTitle("N(event)");             //set yaxis(n event)

    // Y axis histo[0] plot settings
    histo[0]->GetYaxis()->SetTitleSize(18);
    histo[0]->GetYaxis()->SetTitleFont(43);
    histo[0]->GetYaxis()->SetLabelSize(7.5);
    histo[0]->GetYaxis()->SetLabelFont(43);
    histo[0]->GetYaxis()->SetTitleOffset(1.3);
    histo[0]->GetXaxis()->SetRangeUser(xmin[0],xmax[0]);   //set x range
    histo[0]->Draw("e");
    can->Write();
    can->Print((hist[0]+".pdf").c_str());   //save pdf
    
    histo[1]->SetLineColor(2);
    histo[1]->SetLineWidth(2);
    histo[1]->SetName("Ele1_Pt_Ele35");
    histo[1]->Draw("samese"); //histsames + e, draw in same histogram
    
    histo[2]->SetLineColor(3);
    histo[2]->SetLineWidth(2);
    histo[2]->SetName("Ele1_Pt_Ele35/Photon200");
    histo[2]->Draw("samese");
    
    histo[3]->SetLineColor(4);
    histo[3]->SetLineWidth(2);
    histo[3]->SetName("Ele1_Pt_Ele35/Photon200/Ele115");
    histo[3]->Draw("samese");
    
    gPad->Update();
    TPaveStats *st0 = (TPaveStats*)histo[0]->FindObject("stats");
    st0->SetY1NDC(0.79); //new y end position
    st0->SetY2NDC(0.69); //new y start position
    st0->SetX1NDC(0.83); //new x start position
    st0->SetX2NDC(1); //new x end position
    
    gPad->Update();
    TPaveStats *st1 = (TPaveStats*)histo[1]->FindObject("stats");
    st1->SetY1NDC(0.68); //new y end position
    st1->SetY2NDC(0.58); //new y start position
    st1->SetX1NDC(0.83); //new x start position
    st1->SetX2NDC(1); //new x end position
    
    gPad->Update();
    TPaveStats *st2 = (TPaveStats*)histo[2]->FindObject("stats");
    st2->SetY1NDC(0.58); //new y end position
    st2->SetY2NDC(0.48); //new y start position
    st2->SetX1NDC(0.83); //new x start position
    st2->SetX2NDC(1); //new x end position
    
    gPad->Update();
    TPaveStats *st3 = (TPaveStats*)histo[3]->FindObject("stats");
    st3->SetY1NDC(0.48); //new y end position
    st3->SetY2NDC(0.38); //new y start position
    st3->SetX1NDC(0.83); //new x start position
    st3->SetX2NDC(1); //new x end position
    
    //legend settings: new TLegend(x1,y1,x2,y2);
    TLegend* leg1 = new TLegend(0.86,0.95,1,0.79);
    leg1->AddEntry(histo[0],"Ele1_Pt","l");
    leg1->AddEntry(histo[1],"Ele1_Pt_Ele35","l");
    leg1->AddEntry(histo[2],"Ele1_Pt_Ele35/Photon200","l");
    leg1->AddEntry(histo[3],"Ele1_Pt_Ele35/Photon200/Ele115","l");
    leg1->SetBorderSize(0);
    leg1->Draw();


    // lower plot(ratio plot) will be in canvas1-pad2
    can->cd(1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.1);
    pad2->SetBottomMargin(0.1);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    TGraphAsymmErrors *ratioGraph = new TGraphAsymmErrors();
    ratioGraph->Divide(histo[2],histo[3]); //Clopper-Pearson interval

    // Ratio plot (ratioGraph) settings
    ratioGraph->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    ratioGraph->Draw("ap");  //axis point
    ratioGraph->GetYaxis()->SetTitle("no_ele115/with_ele115");
    ratioGraph->GetYaxis()->SetNdivisions(505);
    ratioGraph->GetYaxis()->SetTitleSize(20);
    ratioGraph->GetYaxis()->SetTitleFont(43);
    ratioGraph->GetYaxis()->SetTitleOffset(1);
    ratioGraph->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratioGraph->GetYaxis()->SetLabelSize(11);
    ratioGraph->GetYaxis()->SetRangeUser(0.95,1.02);
    ratioGraph->GetXaxis()->SetLimits(xmin[0],xmax[0]);


    // X axis ratio plot settings
    ratioGraph->GetXaxis()->SetTitleSize(20);
    ratioGraph->GetXaxis()->SetTitleFont(43);
    ratioGraph->GetXaxis()->SetTitleOffset(4.);
    ratioGraph->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratioGraph->GetXaxis()->SetLabelSize(10);
    ratioGraph->Draw("ap");  //axis point

    TLine *l=new TLine(xmin[0],1,xmax[0],1);
    l->SetLineColor(kBlue);
    l->Draw();
    
    can->cd(2);   //go to canvas2
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
    pad3->SetBottomMargin(0.1); // Upper and lower plot are not joined
    pad3->SetGridx();         // Vertical grid
//    pad3->SetLogy();
    pad3->Draw();             // Draw the upper pad: pad1
    pad3->cd();               // pad1 becomes the current pad
    
    histo[4]->SetLineColor(1);
    histo[4]->SetLineWidth(2);
    histo[4]->SetName("Ele2_Pt");         //set name
    histo[4]->GetXaxis()->SetTitle(xtitiles[4].c_str());    //set xaxis
    histo[4]->GetYaxis()->SetTitle("N(event)");             //set yaxis(n event)

    // Y axis histo[0] plot settings
    histo[4]->GetYaxis()->SetTitleSize(18);
    histo[4]->GetYaxis()->SetTitleFont(43);
    histo[4]->GetYaxis()->SetLabelSize(7.5);
    histo[4]->GetYaxis()->SetLabelFont(43);
    histo[4]->GetYaxis()->SetTitleOffset(1.3);
    histo[4]->GetXaxis()->SetRangeUser(xmin[4],xmax[4]);   //set x range
    histo[4]->Draw("e");
    can->Write();
    can->Print((hist[4]+".pdf").c_str());   //save pdf
    
    histo[5]->SetLineColor(2);
    histo[5]->SetLineWidth(2);
    histo[5]->SetName("Ele2_Pt_Ele35");
    histo[5]->Draw("samese"); //histsames + e, draw in same histogram
    
    histo[6]->SetLineColor(3);
    histo[6]->SetLineWidth(2);
    histo[6]->SetName("Ele2_Pt_Ele35/Photon200");
    histo[6]->Draw("samese");
    
    histo[7]->SetLineColor(4);
    histo[7]->SetLineWidth(2);
    histo[7]->SetName("Ele2_Pt_Ele35/Photon200/Ele115");
    histo[7]->Draw("samese");
    
    gPad->Update();
    TPaveStats *st4 = (TPaveStats*)histo[4]->FindObject("stats");
    st4->SetY1NDC(0.79); //new y end position
    st4->SetY2NDC(0.69); //new y start position
    st4->SetX1NDC(0.83); //new x start position
    st4->SetX2NDC(1); //new x end position
    
    gPad->Update();
    TPaveStats *st5 = (TPaveStats*)histo[5]->FindObject("stats");
    st5->SetY1NDC(0.68); //new y end position
    st5->SetY2NDC(0.58); //new y start position
    st5->SetX1NDC(0.83); //new x start position
    st5->SetX2NDC(1); //new x end position
    
    gPad->Update();
    TPaveStats *st6 = (TPaveStats*)histo[6]->FindObject("stats");
    st6->SetY1NDC(0.58); //new y end position
    st6->SetY2NDC(0.48); //new y start position
    st6->SetX1NDC(0.83); //new x start position
    st6->SetX2NDC(1); //new x end position
    
    gPad->Update();
    TPaveStats *st7 = (TPaveStats*)histo[7]->FindObject("stats");
    st7->SetY1NDC(0.48); //new y end position
    st7->SetY2NDC(0.38); //new y start position
    st7->SetX1NDC(0.83); //new x start position
    st7->SetX2NDC(1); //new x end position
    
    //legend settings: new TLegend(x1,y1,x2,y2);
    TLegend* leg2 = new TLegend(0.86,0.95,1,0.79);
    leg2->AddEntry(histo[4],"Ele1_Pt","l");
    leg2->AddEntry(histo[5],"Ele1_Pt_Ele35","l");
    leg2->AddEntry(histo[6],"Ele1_Pt_Ele35/Photon200","l");
    leg2->AddEntry(histo[7],"Ele1_Pt_Ele35/Photon200/Ele115","l");
    leg2->SetBorderSize(0);
    leg2->Draw();


    // lower plot(ratio plot) will be in canvas1-pad2
    can->cd(2);
    TPad *pad4 = new TPad("pad4", "pad4", 0, 0.05, 1, 0.3);
    pad4->SetTopMargin(0.1);
    pad4->SetBottomMargin(0.1);
    pad4->SetGridx(); // vertical grid
    pad4->Draw();
    pad4->cd();       // pad2 becomes the current pad

    TGraphAsymmErrors *ratioGraph2 = new TGraphAsymmErrors();
    ratioGraph2->Divide(histo[6],histo[7]); //print value, if wanna print value and error use pois+v: ratioGraph->Divide(h2,h0,"poisv")

    // Ratio plot (ratioGraph) settings
    ratioGraph2->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    ratioGraph2->Draw("ap");  //axis point
    ratioGraph2->GetYaxis()->SetTitle("without_ele115/with_ele115");
    ratioGraph2->GetYaxis()->SetNdivisions(505);
    ratioGraph2->GetYaxis()->SetTitleSize(20);
    ratioGraph2->GetYaxis()->SetTitleFont(43);
    ratioGraph2->GetYaxis()->SetTitleOffset(1);
    ratioGraph2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratioGraph2->GetYaxis()->SetLabelSize(11);
    ratioGraph2->GetYaxis()->SetRangeUser(0.95,1.02);
    ratioGraph2->GetXaxis()->SetLimits(xmin[4],xmax[4]);


    // X axis ratio plot settings
    ratioGraph2->GetXaxis()->SetTitleSize(20);
    ratioGraph2->GetXaxis()->SetTitleFont(43);
    ratioGraph2->GetXaxis()->SetTitleOffset(4.);
    ratioGraph2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratioGraph2->GetXaxis()->SetLabelSize(10);
    ratioGraph2->Draw("ap");  //axis point

    TLine *l2=new TLine(xmin[4],1,xmax[4],1);
    l2->SetLineColor(kBlue);
    l2->Draw();
    can->Write();
     lqtoxe.Close();
}
