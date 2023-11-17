#include "ParserClass.h"

#include <fstream>
#include <cstdlib>


bool ParserClass::ReadFile(string fname) {
 ifstream infile;
 infile.open(fname.c_str());
 if (!infile.is_open()) return false;
 streampos fileSize = 0;
 infile.seekg(fileSize, ios_base::end);
 fileSize = infile.tellg();
 infile.seekg(0, ios_base::beg);

 Buffer.resize(fileSize);
 infile.read(&(Buffer[0]), fileSize);
 infile.close();
// cerr << "Read " << fileSize << " characters.\n";
 Pos = 0;
 return true;
}


bool ParserClass::FindToken(string token) {
 int toklen = token.size();
 int tokenPos = Buffer.find(token, Pos);
// cerr << "pos = " << tokenPos << endl;
 if (tokenPos == -1)
  return false;
 else {
  Pos = tokenPos+token.size();
  return true;
 }
}


bool ParserClass::ReadInt (int &val) {
// int numChars;
// int success = sscanf (&(Buffer[Pos]), " %d %n", &val, &numChars);
// if (success) Pos += numChars;
// return (success == 1);
 char * endptr;
 val = strtol (&(Buffer[Pos]), &endptr, 10);
 if (endptr == NULL) return false;
 Pos += endptr - &(Buffer[Pos]);
 return true;
}


bool ParserClass::ReadDouble (double &val) {
 char *endptr;
 val =strtod (&(Buffer[Pos]), &endptr);
 if (endptr == NULL) return false;
 Pos += endptr - &(Buffer[Pos]);
 return true;
}


bool ParserClass::ReadComplex (complex<double> &val) {
 double re, im;
 if (FindToken ("("))
  if (ReadDouble(re))
   if (FindToken(","))
    if (ReadDouble(im))
     if (FindToken(")")) {
      val = complex<double>(re,im);
      return true;
     }
 return false;
}


bool isWhiteSpace (char c) {
 return ((c==' ') || (c=='\t') || (c=='\n') || (c=='\r'));
}


bool ParserClass::ReadWord (string &word) {
 bool found = false;
 word = "";
 char str[2];
 str[1] = '\0';
 while (isWhiteSpace (Buffer[Pos]) && (Pos<(Buffer.size()-1))) Pos++;
 while (!isWhiteSpace(Buffer[Pos]) && (Pos<Buffer.size()-1)) {
  str[0] = Buffer[Pos];
  word.append(str);
  found = true;
  Pos++;
 }
 return found;
}
