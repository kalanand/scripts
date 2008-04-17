// sets up the global chains for Qinglin's analysis

#include "../BToDKTo3piK/globals/chains.hh"

void setupChainsQZ() {
  cout << "--- START setupChains() ---" << endl;
// Signal
    chainRanDK.Add( RanfilePathQZ  + "sig-d0k-3pi-1.root");
    chainRanDK.Add( RanfilePathQZ  + "sig-d0k-3pi-2.root");

//Bch sample 
    chainRanBchDpi.Add( RanfilePathQZ  + "bch-d0k-3pi-d0pi-1.root");

    chainRanBchComb.Add( RanfilePathQZ  + "bch-d0k-3pi-other-1.root");
    chainRanBchComb.Add( RanfilePathQZ  + "bch-d0k-3pi-other-2.root");
    chainRanBchComb.Add( RanfilePathQZ  + "bch-d0k-3pi-other-3.root");


//B0 sample
    chainRanB0Comb.Add( RanfilePathQZ  + "b0-d0k-3pi-other-1.root");

//cc sample
    chainRanCc.Add( RanfilePathQZ  + "cc-d0k-3pi-1.root");
    chainRanCc.Add( RanfilePathQZ  + "cc-d0k-3pi-2.root"); 
    chainRanCc.Add( RanfilePathQZ  + "cc-d0k-3pi-3.root"); 

//uds Sample
    chainRanUds.Add( RanfilePathQZ  + "uds-d0k-3pi-1.root");      
    chainRanUds.Add( RanfilePathQZ  + "uds-d0k-3pi-2.root");      
    chainRanUds.Add( RanfilePathQZ  + "uds-d0k-3pi-3.root");      
    chainRanUds.Add( RanfilePathQZ  + "uds-d0k-3pi-4.root");      

// off peak sample 
   chainRanOffData.Add( RanfilePathQZ  + "off-d0k-3pi-1.root");

// on peak Sample
   chainRanOnData.Add( RanfilePathQZ  + "on-d0k-3pi-1.root");
   chainRanOnData.Add( RanfilePathQZ  + "on-d0k-3pi-2.root");
   chainRanOnData.Add( RanfilePathQZ  + "on-d0k-3pi-3.root");
   chainRanOnData.Add( RanfilePathQZ  + "on-d0k-3pi-4.root");
   chainRanOnData.Add( RanfilePathQZ  + "on-d0k-3pi-5.root");
   chainRanOnData.Add( RanfilePathQZ  + "on4-d0k-3pi-1.root");
   chainRanOnData.Add( RanfilePathQZ  + "on4-d0k-3pi-3.root");
   chainRanOnData.Add( RanfilePathQZ  + "on4-d0k-3pi-2.root");
  cout << "--- END setupChains() ---" << endl;
}

