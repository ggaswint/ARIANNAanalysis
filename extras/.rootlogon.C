{
   gSystem->SetIncludePath("-g -DROOTSUPPORT -I$ROOTSYS/include "
			   "-I$SNS/include "
			   "-I./"
			   );

   gSystem->Setenv("TZ","UTC");
   gStyle->SetTimeOffset(0);

   gROOT->SetStyle("Pub");
   //gStyle->SetPalette(1);
   gStyle->SetOptStat(1111111);
   gStyle->SetStatColor(0);
   gStyle->SetLabelColor(1,"xyz");
   gStyle->SetLabelSize(0.05,"xyz");
   gStyle->SetTitleSize(0.05,"xyz");
   gStyle->SetLineWidth(1);
   gStyle->SetHistLineWidth(2);
   gStyle->SetFrameLineWidth(2);
   gStyle->SetFuncWidth(3);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleColor(1);
   gStyle->SetTitleYOffset(1.6);
   gStyle->SetTitleTextColor(1);

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.13);
   gStyle->SetPadTopMargin(0.10);
   gStyle->SetPadBottomMargin(0.15);

   const UInt_t ncolors = 100;
   /*
   // dark blue to azure
   const UInt_t NRGBs = 3;
   Double_t Stops[NRGBs] = { 0.00, 0.60, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.00, 0.00 };
   Double_t Green[NRGBs] = { 0.00, 0.20, 0.90 };
   Double_t Blue[NRGBs]  = { 0.20, 0.90, 1.00 };
   */
   
   // dark blue to red to gold
   const UInt_t NRGBs = 4;
   Double_t Stops[NRGBs] = { 0.00, 0.45, 0.65, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.55, 0.85, 1.00 };
   Double_t Green[NRGBs] = { 0.00, 0.00, 0.00, 0.85 };
   Double_t Blue[NRGBs]  = { 0.25, 0.25, 0.00, 0.00 };
   
   /*
   // dark blue to purple to red
   const UInt_t NRGBs = 5;
   Double_t Stops[NRGBs] = { 0.00, 0.30, 0.65, 0.80, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.00, 0.65, 0.80, 1.00 };
   Double_t Green[NRGBs] = { 0.00, 0.25, 0.10, 0.00, 0.20 };
   Double_t Blue[NRGBs]  = { 0.20, 0.55, 0.55, 0.30, 0.00 };
   */
/*   
   // dark blue to purple to red (a bit brighter)
   const UInt_t NRGBs = 5;
   Double_t Stops[NRGBs] = { 0.00, 0.30, 0.65, 0.80, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.00, 0.75, 0.80, 1.00 };
   Double_t Green[NRGBs] = { 0.00, 0.35, 0.10, 0.00, 0.20 };
   Double_t Blue[NRGBs]  = { 0.25, 0.75, 0.65, 0.30, 0.00 };
*/   
   /*
   // dark blue to purple to red to gold
   const UInt_t NRGBs = 6;
   Double_t Stops[NRGBs] = { 0.00, 0.30, 0.65, 0.80, 0.95, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.00, 0.75, 0.80, 1.00, 1.00 };
   Double_t Green[NRGBs] = { 0.00, 0.35, 0.10, 0.00, 0.20, 0.85 };
   Double_t Blue[NRGBs]  = { 0.25, 0.75, 0.65, 0.30, 0.00, 0.00 };
   */
   /*
   // dark red to gold
   const UInt_t NRGBs = 3;
   Double_t Stops[NRGBs] = { 0.00, 0.55, 1.00 };
   Double_t Red[NRGBs]   = { 0.20, 0.55, 1.00 };
   Double_t Green[NRGBs] = { 0.00, 0.00, 0.85 };
   Double_t Blue[NRGBs]  = { 0.00, 0.00, 0.00 };
   */
   /*
   const UInt_t NRGBs = 3;
   Double_t Stops[NRGBs] = { 0.00, 0.55, 1.00 };
   Double_t Red[NRGBs]   = { 0.30, 0.55, 0.00 };
   Double_t Green[NRGBs] = { 0.00, 0.00, 0.80 };
   Double_t Blue[NRGBs]  = { 0.00, 0.65, 0.90 };
   */
   /*
   // dark blue to teal
   const UInt_t NRGBs = 3;
   Double_t Stops[NRGBs] = { 0.00, 0.50, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.00, 0.20 };
   Double_t Green[NRGBs] = { 0.00, 0.90, 1.00 };
   Double_t Blue[NRGBs]  = { 0.20, 1.00, 0.50 };
   */
   /*
   // turquoise
   const UInt_t NRGBs = 3;
   Double_t Stops[NRGBs] = { 0.00, 0.60, 1.00 };
   Double_t Red[NRGBs]   = { 0.00, 0.00, 0.00 };
   Double_t Green[NRGBs] = { 0.00, 0.60, 0.95 };
   Double_t Blue[NRGBs]  = { 0.10, 0.60, 0.80 };
   */
   /*
   // rainbow with from deep purple to deep red
   const UInt_t NRGBs = 8;
   Double_t Stops[NRGBs] = { 0.00, 0.07, 0.15, 0.35, 0.50, 0.70, 0.90, 1.00 };
   Double_t Red[NRGBs]   = { 0.20, 0.30, 0.00, 0.00, 0.00, 0.97, 0.97, 0.60 };
   Double_t Green[NRGBs] = { 0.00, 0.04, 0.00, 0.97, 0.97, 0.97, 0.00, 0.00 };
   Double_t Blue[NRGBs]  = { 0.20, 0.30, 0.50, 0.97, 0.00, 0.00, 0.00, 0.00 };
   */
#if ROOT_VERSION_CODE < ROOT_VERSION(5,15,4)
   TStyle::CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, ncolors);
#else
   TColor::InitializeColors();
   TColor::CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, ncolors);
#endif
}
