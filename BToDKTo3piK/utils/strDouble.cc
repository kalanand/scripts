
// Turns a double into a string without the stupid spaces

TString strDouble(double d) {
  TString num;
  num += d;
  
  const char * nameStr = num;
  while(nameStr[0] == ' ') {
    ++nameStr;
  }
  
  TString result = nameStr;
  return result;
}
  
