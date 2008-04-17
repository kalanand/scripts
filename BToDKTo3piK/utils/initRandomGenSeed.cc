#include <fstream.h>
#include "../BToDKTo3piK/utils/setRandomGenSeed.cc"

UInt_t initRandomGenSeed(const char* seedFile) {

   UInt_t seed(63287);
   Int_t seedValue(0);
   Bool_t found = kFALSE;
   ifstream ifs(seedFile);  
 
   RooStreamParser parser(ifs);
   parser.setPunctuation("=");
   TString token;

   while ( !(ifs.eof() || ifs.fail() || parser.atEOF()) ) {

        // Read next token until end of file
        token = parser.readToken() ;

        // Skip empty or comment lines
        if ( token.IsNull() || token.Contains("//") ) {
            continue ;
        }

        // check whether we have the correct token here
        if ( token == "seedValue" && !parser.expectToken("=",kTRUE) ) {
            token = parser.readToken();
            if ( ! token.IsNull() ) {
                if ( !parser.convertToInteger(token,seedValue) ) {
                    seed = seedValue;
                    found = kTRUE;
                }
            }
        }
        parser.zapToEnd();
    }
    ifs.close();

    // set the new random seed
    setRandomGenSeed(seed);

    // increment seed value and write to file
    seed++;
    ofstream ofs(seedFile);
    if ( ofs.fail() ) {
        cout << "FitUtils::initRandomGenSeed: "
             << " cannot write to " << seedFile << endl;
    } else {
        ofs << "//" << endl;
        ofs << "// Random generator seed file" << endl;
        ofs << "// This file is updated automagically." << endl;
        ofs << "// Do not edit!" << endl;
        ofs << "//" << endl;
        ofs << "   seedValue = " << seed << endl;
        ofs << endl;
        ofs.close();
    }
    return --seed;
}
