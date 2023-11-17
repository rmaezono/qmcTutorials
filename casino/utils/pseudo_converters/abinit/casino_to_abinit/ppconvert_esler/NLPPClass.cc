#include "NLPPClass.h"
#include "XMLWriterClass2.h"
#include "ParserClass.h"
#include "ParseCommand.h"
#include <ctime>
#include <sstream>


void ChannelPotentialClass::WriteChannel (XMLWriterClass &writer,
 bool writeVl) {
 string channels[] = {"s", "p", "d", "f", "g", "h", "i", "j"};
 // Use a logarithmic grid:
 // r(i) = a*(exp(b*i) - 1)
 const int numPoints = 2002;
 double end = Vl.grid.End();
 double step  = 0.00625;
 double scale = end/(exp(step*(numPoints-1))-1.0);
 if (writeVl)
  writer.StartElement("vps");
 else
  writer.StartElement("pswf");
 writer.WriteAttribute("principal-n", n_principal);
 writer.WriteAttribute("l", channels[l]);
 writer.WriteAttribute("spin", -1);
 if (writeVl) {
  writer.WriteAttribute("cutoff", Cutoff);
  writer.WriteAttribute("occupation", Occupation);
 }
 writer.StartElement("radfunc");
 writer.StartElement("grid");
 writer.WriteAttribute("type", "log");
 writer.WriteAttribute("units", "bohr");
 writer.WriteAttribute("scale", scale, true);
 writer.WriteAttribute("step", step, true);
 writer.WriteAttribute("npts", numPoints-1);
 writer.EndElement(); // "grid"
 vector<double> data;
 for (int i=1; i<numPoints; i++) {
  double r  = scale * (exp(step*i)-1.0);
  if (writeVl) {
   r = min (r, Vl.grid.End());
   data.push_back(Vl(r));
  } else {
   r = min (r, ul.grid.End());
   data.push_back(ul(r));
  }
 }
 writer.WriteElement("data", data);
 writer.EndElement(); // "radfunc"
 writer.EndElement(); // "vps" or "pswf"
}


void NLPPClass::SetupMaps() {
 UnitToHartreeMap[string("hartree")] = 1.0;
 UnitToHartreeMap[string("rydberg")] = 0.5;
 UnitToHartreeMap[string("ev")]      = 0.03674932595264097934;

 UnitToBohrMap[string("bohr")]     = 1.0;
 UnitToBohrMap[string("atomic")]   = 1.0;
 UnitToBohrMap[string("angstrom")] = 1.8897261;

 XCMap[XC_LDA] ="LDA" ; XCRevMap["LDA"] =XC_LDA;  XCRevMap["lda"] =XC_LDA;
 XCMap[XC_GGA] ="GGA" ; XCRevMap["GGA"] =XC_GGA;  XCRevMap["gga"] =XC_LDA;
 XCMap[XC_HF]  ="HF"  ; XCRevMap["HF"]  =XC_HF;   XCRevMap["hf"]  =XC_LDA;
 XCMap[XC_DF]  ="DF"  ; XCRevMap["DF"]  =XC_DF;   XCRevMap["df"]  =XC_DF;
 XCMap[XC_NONE]="NONE"; XCRevMap["NONE"]=XC_NONE; XCRevMap["none"]=XC_NONE;

 ChannelMap[0] = "s";  ChannelRevMap["s"] = 0;
 ChannelMap[1] = "p";  ChannelRevMap["p"] = 1;
 ChannelMap[2] = "d";  ChannelRevMap["d"] = 2;
 ChannelMap[3] = "f";  ChannelRevMap["f"] = 3;
 ChannelMap[4] = "g";  ChannelRevMap["g"] = 4;
 ChannelMap[5] = "h";  ChannelRevMap["h"] = 5;
 ChannelMap[6] = "i";  ChannelRevMap["i"] = 6;
 ZToSymbolMap[1]   = "H";  ZToSymbolMap[2]   = "He";
 ZToSymbolMap[3]   = "Li"; ZToSymbolMap[4]   = "Be";
 ZToSymbolMap[5]   = "B";  ZToSymbolMap[6]   = "C";
 ZToSymbolMap[7]   = "N";  ZToSymbolMap[8]   = "O";
 ZToSymbolMap[9]   = "F";  ZToSymbolMap[10]  = "Ne";
 ZToSymbolMap[11]  = "Na"; ZToSymbolMap[12]  = "Mg";
 ZToSymbolMap[13]  = "Al"; ZToSymbolMap[14]  = "Si";
 ZToSymbolMap[15]  = "P";  ZToSymbolMap[16]  = "S";
 ZToSymbolMap[17]  = "Cl"; ZToSymbolMap[18]  = "Ar";
 ZToSymbolMap[19]  = "K";  ZToSymbolMap[20]  = "Ca";
 ZToSymbolMap[21]  = "Sc"; ZToSymbolMap[22]  = "Ti";
 ZToSymbolMap[23]  = "V";  ZToSymbolMap[24]  = "Cr";
 ZToSymbolMap[25]  = "Mn"; ZToSymbolMap[26]  = "Fe";
 ZToSymbolMap[27]  = "Co"; ZToSymbolMap[28]  = "Ni";
 ZToSymbolMap[29]  = "Cu"; ZToSymbolMap[30]  = "Zn";
 ZToSymbolMap[31]  = "Ga"; ZToSymbolMap[32]  = "Ge";
 ZToSymbolMap[33]  = "As"; ZToSymbolMap[34]  = "Se";
 ZToSymbolMap[35]  = "Br"; ZToSymbolMap[36]  = "Kr";
 ZToSymbolMap[37]  = "Rb"; ZToSymbolMap[38]  = "Sr";
 ZToSymbolMap[39]  = "Y";  ZToSymbolMap[40]  = "Zr";
 ZToSymbolMap[41]  = "Nb"; ZToSymbolMap[42]  = "Mo";
 ZToSymbolMap[43]  = "Tc"; ZToSymbolMap[44]  = "Ru";
 ZToSymbolMap[45]  = "Rh"; ZToSymbolMap[46]  = "Pd";
 ZToSymbolMap[47]  = "Ag"; ZToSymbolMap[48]  = "Cd";
 ZToSymbolMap[49]  = "In"; ZToSymbolMap[50]  = "Sn";
 ZToSymbolMap[51]  = "Sb"; ZToSymbolMap[52]  = "Te";
 ZToSymbolMap[53]  = "I";  ZToSymbolMap[54]  = "Xe";
 ZToSymbolMap[55]  = "Cs"; ZToSymbolMap[56]  = "Ba";
 ZToSymbolMap[57]  = "La"; ZToSymbolMap[58]  = "Ce";
 ZToSymbolMap[59]  = "Pr"; ZToSymbolMap[60]  = "Nd";
 ZToSymbolMap[61]  = "Pm"; ZToSymbolMap[62]  = "Sm";
 ZToSymbolMap[63]  = "Eu"; ZToSymbolMap[64]  = "Gd";
 ZToSymbolMap[65]  = "Tb"; ZToSymbolMap[66]  = "Dy";
 ZToSymbolMap[67]  = "Ho"; ZToSymbolMap[68]  = "Er";
 ZToSymbolMap[69]  = "Tm"; ZToSymbolMap[70]  = "Yb";
 ZToSymbolMap[71]  = "Lu"; ZToSymbolMap[72]  = "Hf";
 ZToSymbolMap[73]  = "Ta"; ZToSymbolMap[74]  = "W";
 ZToSymbolMap[75]  = "Re"; ZToSymbolMap[76]  = "Os";
 ZToSymbolMap[77]  = "Ir"; ZToSymbolMap[78]  = "Pt";
 ZToSymbolMap[79]  = "Au"; ZToSymbolMap[80]  = "Hg";
 ZToSymbolMap[81]  = "Tl"; ZToSymbolMap[82]  = "Pb";
 ZToSymbolMap[83]  = "Bi"; ZToSymbolMap[84]  = "Po";
 ZToSymbolMap[85]  = "At"; ZToSymbolMap[86]  = "Rn";
 ZToSymbolMap[87]  = "Fr"; ZToSymbolMap[88]  = "Ra";
 ZToSymbolMap[89]  = "Ac"; ZToSymbolMap[90]  = "Th";
 ZToSymbolMap[91]  = "Pa"; ZToSymbolMap[92]  = "U";
 ZToSymbolMap[93]  = "Np"; ZToSymbolMap[94]  = "Pu";
 ZToSymbolMap[95]  = "Am"; ZToSymbolMap[96]  = "Cm";
 ZToSymbolMap[97]  = "Bk"; ZToSymbolMap[98]  = "Cf";
 ZToSymbolMap[99]  = "Es"; ZToSymbolMap[100] = "Fm";
 ZToSymbolMap[101] = "Mc"; ZToSymbolMap[102] = "No";
 ZToSymbolMap[103] = "Lw";
}


void NLPPClass::WriteXML(string fname) {
// int rc;
// xmlTextWriterPtr writer = xmlNewTextWriterFilename (fname.c_str(), 0);

// rc = xmlTextWriterStartDocument (writer, NULL, "UTF-8", NULL);
// // Start pseudo
// rc = xmlTextWriterStartElement(writer, (xmlChar*)"pseudo");
// rc = xmlTextWriterWriteAttribute (writer, (xmlChar*) "version",
//  (xmlChar*)"0.5");
// rc = xmlTextWriterEndElement(writer); // "pseudo"
// rc = xmlTextWriterEndDocument (writer);
// xmlFreeTextWriter(writer);
 for (int l=0; l<ChannelPotentials.size(); l++)
  if (!ChannelPotentials[l].HasProjector) {
   cerr << "Missing projector for l=" << l << ".  Aborting." << endl;
   exit(1);
  }
 XMLWriterClass writer;
 writer.StartDocument(fname, "1.0", "UTF-8");
 writer.StartElement("pseudo");
 writer.WriteAttribute ("version", "0.5");
 writer.StartElement("header");
 writer.WriteAttribute("symbol", ZToSymbolMap[AtomicNumber]);
 writer.WriteAttribute("atomic-number", (double)AtomicNumber);
 writer.WriteAttribute("zval", PseudoCharge);
 writer.WriteAttribute("relativistic", Relativistic ? "yes" : "no");
 writer.WriteAttribute("polarized", "no");
 writer.WriteAttribute("creator", "ppconvert");
 writer.WriteAttribute("flavor", "Troullier-Martins");

 writer.WriteAttribute("core-corrections", "no");
// writer.WriteAttribute("xc-functional-type", XCMap[XC]);
 writer.WriteAttribute("xc-functional-type", "GGA");
 writer.WriteAttribute("xc-functional-parametrization",
  "Perdew-Burke-Ernzerhof");
 writer.EndElement(); // "header"

 // Write the grid information:
 const int numPoints = 2002;
 double end = PotentialGrid.End();
 double step  = 0.00625;
 double scale = end/(exp(step*(numPoints-1))-1.0);
 writer.StartElement("grid");
 writer.WriteAttribute("type", "log");
 writer.WriteAttribute("units", "bohr");
 writer.WriteAttribute("scale", scale, true);
 writer.WriteAttribute("step", step, true);
 writer.WriteAttribute("npts", numPoints-1);
 writer.EndElement(); // "grid"

 writer.StartElement("semilocal");
 writer.WriteAttribute("units", "hartree");
 writer.WriteAttribute("format", "r*V");
 writer.WriteAttribute("npots-down", (int)ChannelPotentials.size());
 writer.WriteAttribute("npots-up"  , 0);
 writer.WriteAttribute("l-local", LocalChannel);
 for (int l=0; l<ChannelPotentials.size(); l++)
  ChannelPotentials[l].WriteChannel (writer, true);
 writer.EndElement(); // "semilocal"

 writer.StartElement("pseudowave-functions");
 writer.WriteAttribute("units", "electrons/bohr^(-3/2)");
 writer.WriteAttribute("format", "u_n,l (r) = 1/r R_n,l (r)");
 writer.WriteAttribute("n-pseudowave-functions-down",
  (int)ChannelPotentials.size());
 writer.WriteAttribute("n-pseudowave-functions-up", 0);
 for (int l=0; l<ChannelPotentials.size(); l++)
  ChannelPotentials[l].WriteChannel (writer, false);
 writer.EndElement(); // "pseudowave-functions"

 writer.EndElement(); // "pseudo"
 writer.EndDocument();
}


void NLPPClass::WriteABINIT (string fileName) {
 if (fileName == "")
  fileName = ZToSymbolMap[AtomicNumber] + ".psp";

 FILE *fout = fopen (fileName.c_str(), "w");
 assert (fout != NULL);

 time_t now = time(NULL);
 struct tm &ltime = *(localtime (&now));
 stringstream date;
 int year = ltime.tm_year % 100;
 int mon  = ltime.tm_mon + 1;
 int day  = ltime.tm_mday;
 date << ((year<10) ? "0" : "") << year;
 date << ((mon<10)  ? "0" : "") << mon;
 date << ((day<10)  ? "0" : "")  << day;

 fprintf (fout, "%s  %s", ZToSymbolMap[AtomicNumber].c_str(),
  asctime(&ltime));
 fprintf (fout, "%9.5f %9.5f %s"
  "                    zatom, zion, pspdat\n", (double)AtomicNumber,
  PseudoCharge, date.str().c_str());
 const int pspcod = 1;
 const int pspxc  = 0;
 const int lmax   = ChannelPotentials.size()-1;
 // HACK HACK HAK
 const int lloc   = LocalChannel;
 const int mmax  = 2001;
 double r2well = 0.0;
 // The following are informational only
 const double e99    = 0.0;
 const double e999   = 0.0;
 const double rms    = 0.0;
 const double ekb1   = 0.0;
 const double ekb2   = 0.0;
 const double epsatm = 0.0;
 const double fchrg  = 0.0;
 const double rchrg  = 1.0;
 const double qchrg  = 0.0;

 fprintf (fout, "%d   %d   %d   %d   %d   %8.5f"
  "              pspcod,pspxc,lmax,lloc,mmax,r2well\n",
  pspcod, pspxc, lmax, lloc, mmax, r2well);

 for (int l=0; l<ChannelPotentials.size(); l++) {
  // The local channel needs no projector
  int nproj = (l == lloc) ? 0 : 1;
  fprintf (fout, "%d %5.8f %5.8f %d %5.8f"
   "          l,e99.0,e99.9,nproj,rcpsp\n",
   l, e99, e999, nproj, ChannelPotentials[l].Cutoff);
  fprintf (fout, "%14.10f %14.10f %14.10f %14.10f"
   "  rms,ekb1,ekb2,epsatm\n", rms, ekb1, ekb2, epsatm);
 }
 fprintf (fout, "%8.5f %8.5f %8.5f"
  "                    rchrg,fchrg,qchrg\n", rchrg, fchrg, qchrg);

 // Write out potentials
 for (int l=0; l<ChannelPotentials.size(); l++) {
  fprintf (fout, "%d = l for CASINO pseudopotential\n", l);
  for (int i=0; i<mmax; i++) {
   double x = (double)i/(double)(mmax-1);
   x += 0.01;
   double r = 100.0*x*x*x*x*x - 1.0e-8;
   double Vl = ChannelPotentials[l].Vl(r)/r;
   fprintf (fout, "%23.16e ", Vl);
   if ((i%3)==2) fprintf (fout, "\n");
  }
 }
 // Now write out radial wave functions
 for (int l=0; l<ChannelPotentials.size(); l++) {
  fprintf (fout, "%d = l for CASINO wave function\n", l);
  for (int i=0; i<mmax; i++) {
   double x = (double)i/(double)(mmax-1);
   x += 0.01;
   double r = 100.0*x*x*x*x*x - 1.0e-8;
   double ul = ChannelPotentials[l].ul(r);
   fprintf (fout, "%23.16e ", ul);
   if ((i%3)==2) fprintf (fout, "\n");
  }
 }
 fclose (fout);
}


void NLPPClass::WriteASCII() {
 FILE *fout = fopen ("pp.dat", "w");
 for (int i=1; i<PotentialGrid.NumPoints(); i++) {
  double r= PotentialGrid[i];
  fprintf (fout, "%24.16e ", r);
  for (int l=0; l<ChannelPotentials.size(); l++)
   fprintf (fout, "%24.16e ", ChannelPotentials[l].Vl(r)/r);
  fprintf (fout, "\n");
 }
 fclose (fout);
}


void NLPPClass::ReadCASINO_PP (string fileName) {
 ParserClass parser;
 parser.ReadFile (fileName);

 assert (parser.FindToken("pseudo-charge"));
 assert (parser.ReadInt (AtomicNumber));
 assert (parser.ReadDouble (PseudoCharge));
 assert (parser.FindToken("(rydberg/hartree/ev):"));
 assert (parser.ReadWord(EnergyUnit));
 cerr << "EnergyUnit = " << EnergyUnit << endl;
 assert (parser.FindToken("(0=s,1=p,2=d..)"));
 assert (parser.ReadInt (LocalChannel));
 assert (parser.FindToken("grid points"));
 int numGridPoints;
 assert (parser.ReadInt(numGridPoints));
 vector<double> gridPoints(numGridPoints), Vl(numGridPoints);
 assert (parser.FindToken ("in"));
 assert (parser.ReadWord(LengthUnit));
 cerr << "LengthUnit = " << LengthUnit << endl;
 assert (parser.FindToken("\n"));
 // assert (parser.FindToken("atomic units"));
 for (int i=0; i<numGridPoints; i++) {
  assert (parser.ReadDouble(gridPoints[i]));
  gridPoints[i] *= UnitToBohrMap [LengthUnit];
 }
 PotentialGrid.Init (gridPoints);
 int l=0;
 bool done(false);
 while (!done) {
  if (!parser.FindToken("(L="))
   done = true;
  else {
   assert (parser.ReadInt(l));
   assert (parser.FindToken("\n"));
   for (int i=0; i<numGridPoints; i++) {
    assert (parser.ReadDouble(Vl[i]));
    Vl[i] *= UnitToHartreeMap[EnergyUnit];
   }
   ChannelPotentials.push_back(ChannelPotentialClass());
   ChannelPotentials[l].Vl.Init(PotentialGrid, Vl);
   ChannelPotentials[l].l = l;
  }
 }
 // Now read the summary file to get core radii
 if (!parser.ReadFile ("summary.txt")) {
  cerr << "Could not find summary.txt file.  Aborting.\n";
  abort();
 }
 assert (parser.FindToken("core radii"));
 assert (parser.FindToken("\n"));
 assert (parser.FindToken("\n"));
 // Read core radii
 double rc;
 for (int i=0; i<=l; i++) {
  assert (parser.FindToken (ChannelMap[i]));
  assert (parser.ReadDouble (rc));
  ChannelPotentials[i].Cutoff = rc;
 }

 assert (parser.FindToken("r_loc"));
 assert (parser.FindToken("\n"));
 assert (parser.FindToken("\n"));
 assert (parser.FindToken("(Grid)"));
 double rloc;
 assert (parser.ReadDouble(rloc));
// for (int i=0; i<=l; i++)
//  ChannelPotentials[i].Cutoff = rloc;

 cerr << "Found " << (l+1) << " l-channel potentials.\n";
}


/// Note:  the pseudopotentials must be read before the wave
/// functions.
void NLPPClass::ReadCASINO_WF (string fileName, int l) {
 // Make sure that l is not too high
 assert (l < ChannelPotentials.size());
 ParserClass parser;
 parser.ReadFile (fileName);

 // Make sure this a wave function file
 assert (parser.FindToken ("wave function"));
 // Find atomic number
 assert (parser.FindToken ("Atomic number"));
 int atomicNumber;
 assert (parser.ReadInt(atomicNumber));
 assert (atomicNumber == AtomicNumber);
 int numOrbitals;
 assert (parser.FindToken("number of orbitals"));
 assert (parser.ReadInt(numOrbitals));
 assert (parser.FindToken("Radial grid"));
 assert (parser.FindToken("\n"));
 int numPoints;
 assert (parser.ReadInt(numPoints));
 assert (PotentialGrid.NumPoints() == numPoints);
 // Make sure we have the same grid as the potential file
 for (int i=0; i<numPoints; i++) {
  double r;
  assert(parser.ReadDouble(r));
  assert (fabs(r-PotentialGrid[i]) < 1.0e-10);
 }
 bool orbFound = false;
 for (int i=0; i<numOrbitals; i++) {
  assert (parser.FindToken ("Orbital #"));
  assert (parser.FindToken ("\n"));
  int spin, n, thisl;
  assert (parser.ReadInt(spin));
  assert (parser.ReadInt(n));
  assert (parser.ReadInt(thisl));
  if (thisl == l) {
   vector<double> ul(numPoints);
   for (int i=0; i<numPoints; i++) assert (parser.ReadDouble(ul[i]));
   ChannelPotentials[l].ul.Init(PotentialGrid, ul);
   ChannelPotentials[l].HasProjector = true;
   ChannelPotentials[l].n_principal = n;
   orbFound = true;
  }
 }
 if (!orbFound) {
  cerr << "Could not file orbital with l=" << l
   << " in file ""filename""" << endl;
  exit(1);
 }
}


int NLPPClass::GetNumChannels() { return ChannelPotentials.size(); }


main(int argc, char **argv) {
 // Create list of acceptable command-line parameters
 list<ParamClass> argList;
 argList.push_back(ParamClass("casino_pot", true));
 argList.push_back(ParamClass("casino_us", true));
 argList.push_back(ParamClass("casino_up", true));
 argList.push_back(ParamClass("casino_ud", true));
 argList.push_back(ParamClass("casino_uf", true));
 argList.push_back(ParamClass("xml", true));
 argList.push_back(ParamClass("tm", true));

 CommandLineParserClass parser(argList);
 bool success = parser.Parse(argc, argv);
 if (!success || parser.NumFiles()!=0 || !parser.Found("casino_pot")) {
  cerr << "Usage:  ppconvert --casino_pot=fname [--casino_us=fname "
   << "[--casino_up=fname...]] [--xml fname.xml] [--tm fname.tm]\n";
  exit(1);
 }
 string xmlFile, tmFile;
 string potFile = parser.GetArg("casino_pot");

 NLPPClass nlpp;
 nlpp.ReadCASINO_PP(potFile);
 int numChannels = nlpp.GetNumChannels();
 if (numChannels > 0) {
  assert (parser.Found("casino_us"));
  nlpp.ReadCASINO_WF(parser.GetArg("casino_us"), 0);
 }
 if (numChannels > 1) {
  assert (parser.Found("casino_up"));
  nlpp.ReadCASINO_WF(parser.GetArg("casino_up"), 1);
 }
 if (numChannels > 2) {
  assert (parser.Found("casino_ud"));
  nlpp.ReadCASINO_WF(parser.GetArg("casino_ud"), 2);
 }
 if (numChannels > 3) {
  assert (parser.Found("casino_uf"));
  nlpp.ReadCASINO_WF(parser.GetArg("casino_uf"), 3);
 }

 if (parser.Found("xml")) nlpp.WriteXML(parser.GetArg("xml"));
 if (parser.Found("tm")) nlpp.WriteABINIT(parser.GetArg("tm"));

// nlpp.WriteXML(xmlFile);
// nlpp.WriteABINIT();
 nlpp.WriteASCII();
// nlpp.ReadCASINO_PP ("b_pp.data.HF");
// nlpp.ReadCASINO_WF ("awfn.data_s2p1_2P", 0);
// nlpp.ReadCASINO_WF ("awfn.data_s2p1_2P", 1);
// nlpp.ReadCASINO_WF ("awfn.data_s2d1_2D", 2);
// nlpp.WriteXML("test.xml");
}
