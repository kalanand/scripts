      program FindBestCand 

      IMPLICIT NONE

      include "newmc.inc"

      REAL HMEMOR
      INTEGER HSPACE,IQUEST
      PARAMETER (HSPACE = 25000000)
      COMMON/PAWC/HMEMOR(HSPACE)
      integer LREC,IERR,IOFF,ICYCLE,ISTAT,NCAND
      integer kpimode, trusig, dumpfile
      integer ncat, ncatt 
      integer i, ii, jj, mnum  
      integer idin, idout 
      integer ncomb, ncomb1, ncomb2, ncomb3, ncomb4

      real enkch, thckch, masskch
      real mdpdg, mbpdg, masspi
      real smd1, smd2  ! sigma of D0 for 4 channels
      real smes ! sigma of Mes for 4 channels
      real calchxsq,bestchisq
          
      character TAG1*1500, TAG2*1500, TAG3*1500,TAG4*1500
      parameter(mnum=2)
      parameter(masskch=0.49368)
      parameter(masspi=0.13957)

c smd1 Ks mode, smd2 pi0 mode

      parameter(smd1=0.0078,smd2=0.0133)
      parameter(smes=0.0026)
      parameter(mdpdg=1.8645,mbpdg=5.279)

      integer CanIndex
      integer index(150), rndVal
      
      
      character*110 outfilenm, infilename
      character*110 infilelist,fname
      character*35 pref

c     file indetifier parameters

      integer nwrite, iwrt, nfile
      integer np1, np2
      data iwrt/0/, np1/0/, np2/0/      

c neural net variables:
      real  hlay(25), hlaytmp
      data  hlay/25*0.0/, hlaytmp/0.0/

      real hlaybk(25)
      data hlaybk/25*0./
      real act(35),   bias(35)
      real actbk(35), biasbk(35)
      data act/35*.0/, bias/35*.0/, actbk/35*.0/, biasbk/35*.0/
      
c ********** coefficients for pi0 modes   **********
c suppress continuum
c input level
      real  hlaycoef(10,25) ! read from a txt file coef_i_kkpi0.txt
c out put level
      real  olaycoef(25)      ! read from a txt file coef_o_kkpi0.txt
      data hlaycoef/250*0./, olaycoef/25*0./
   
c suppress all bk
      real hlaybkcoef(9,24) 
      real olaybkcoef(24)
      data hlaybkcoef/216*0./, olaybkcoef/24*0./

      integer tmpnum
c 
c event block
c
      integer myeventNumber, myrunnumber, myUpper
      integer myB1decMode, myB2decMode, myLower   
      integer myD1CaDecMode, myD2CaDecMode 
 
      real myR2, myHem1Mass, myHem2Mass  

      common /tag1/myeventNumber,myrunNumber,myR2,
     +       myLower, myUpper,myB1decMode,
     +       myB2decMode,myD1CaDecMode,myD2CaDecMode,
     +       myHem1Mass, myHem2Mass, ncomb   

      DATA 
     +     myeventNumber/-100/, myrunNumber/-100/,
     +     myR2/-100./,myHem1Mass/-100.0/,
     +     myHem2Mass/-100.0/       

      DATA tag1/'
     +     eventNumber:i,runNumber:i,R2:r,
     +     Lower:i,Upper:i,B1decMode:i, B2decMode:i,
     +     D1CaDecMode:i, D2CaDecMode:i,
     +     Hem1Mass:r,Hem2Mass:r, ncomb:i'/
	
c
cD0Kch block
c
      integer  myKbLepMassInt, myprobB0Lep 
      integer  mydksVtxVctChisqInt, myKsDecLenInt

      integer inp, inhl, inkkp, olkkp, in3pi, ol3pi, ab 
      integer inbkkkp, olbkkkp
      real nntest

      real  
c Bachelor kaon
     +     myHdTrkThc,myHdTrkThcErr,myHdTrkPcm,
     +     myHdTrkKaonNN,myHdTrkCMTheta,
c Ks 
     +     myKsMass,myKsVtxFitChisq,myKsDecLen,mycosDKsVtx,
     +     myKsHel, mydksVtxVctChisq,
c pi0
     +     mypi0Mass,mypi0EAsy,mypi0Hel,mypi0CMmom,
     +     myBest2gmaMass1,myBest2gmaMass2,myBestcos2gmaHel1,
     +     myBestcos2gmaHel2,myBest2gPcm1,myBest2gPcm2, 
c D0 
     +     myd0kkMCMass,myd0piKpMCMass,myd0piKmMCMass,
     +     myd0ppMCMass, myd0pPpMCMass,myd0pPmMCMass,
     +     myd0kkMass,myd0piKpMass,myd0piKmMass,
     +     myd0ppMass, myd0pPpMass,myd0pPmMass,
     +     myd0kkUPMass,myd0piKpUPMass,myd0piKmUPMass,
     +     myd0ppUPMass, myd0pPpUPMass,myd0pPmUPMass, 
     +     myd0VtxFitChisq,
     +     myd0Mass,myd0FkPromptMass,myd0Pcm,myd0fakeMasspKp,
     +     myd0fakeMasspKm,myd0swpmass,myd0flightDist,
c B variables
     +     myBFCLeo,myBFLegendre, myCosBThetaT, 
     +     myKbKUpmass, myKbKlowMass, myKbLepMass,
     +     myKbKsUpMass, myKbKslowMass,myQHemiDiff, myKroeBvtxDoca,
     +     mymixNeuMass, mymixChgMass,myKbD0Doca,myKbD0DocaInlog,
     +     mycosBmomDThrDFm,myBCMcosTheta,mybVtxFitChisq,
     +     mycosBDVtxDmom, myBandBRoeZ,myBandBRoeZInlog, myBandBRoeY,
     +     myprobB0E,myprobB0Mu, mybdVtxVctChisq,     
c beam
     +     mymES,  myDeltaE,
c neuron 
     +     mynnout, mybknnout

      integer myLength, 
c kaon
     +     myHdTrkChge,myHdtrktruth,
c ks
     +     myorigKs,
c pi0
     +     myorigPi0,
c D0
     +     myd0trueKsPP,
     +     myd0DecMode,myKpPidBit,myKmPidBit,myPipPidBit,
     +     myPimPidBit,myTrueD0flg,myrecD0CaDec,myd0RecDec,
c B
     +     myexclTruth,mytrueBflg,myDecInChrm


      common /tag2/ myLength,
     + myHdTrkThc(mnum),         myHdTrkThcErr(mnum), 
     + myHdTrkKaonNN(mnum),      myHdTrkCMTheta(mnum),
     + myHdTrkPcm(mnum),         myHdTrkChge(mnum),
     + myHdTrkTruth(mnum),
     + myKsMass(mnum),           myKsVtxFitChisq(mnum),         
     + myKsDecLen(mnum),         myKsDecLenInt(mnum),   
     + mycosDKsVtx(mnum),
     + myKsHel(mnum),            mydksVtxVctChisq(mnum),
     + mydksVtxVctChisqInt(mnum),myorigKs(mnum),
     + mypi0Mass(mnum),          mypi0EAsy(mnum),
     + mypi0Hel(mnum),           mypi0CMmom(mnum),
     + myBest2gmaMass1(mnum),    myBest2gmaMass2(mnum),
     + myBestcos2gmaHel1(mnum),  myBestcos2gmaHel2(mnum),
     + myBest2gPcm1(mnum),       myBest2gPcm2(mnum),
     + myOrigPi0(mnum)
 
c     DATA myHdTrkThcErr/mnum*(-1000.0)/	

      DATA tag2/' Length[0,2]:i,
     + HdTrkThc(Length):r,         HdTrkThcErr(Length):r,
     + HdTrkKaonNN(Length):r,      HdTrkCMTheta(Length):r,
     + HdTrkPcm(Length):r,         HdTrkChge(Length):i,
     + HdTrkTruth(Length):i,
     + KsMass(Length):r,           KsVtxFitChisq(Length):r,         
     + KsDecLen(Length):r,         KsDecLenInt(Length):i,
     + cosDKsVtx(Length):r,
     + KsHel(Length):r,            dksVtxVctChisq(Length):r,
     + dksVtxVctChisqInt(Length):i,origKs(Length):i,
     + pi0Mass(Length):r,          pi0EAsy(Length):r,
     + pi0Hel(Length):r,           pi0CMmom(Length):r,
     + Best2gmaMass1(Length):r,    Best2gmaMass2(Length):r,
     + Bestcos2gmaHel1(Length):r,  Bestcos2gmaHel2(Length):r,
     + Best2gPcm1(Length):r,       Best2gPcm2(Length):r,
     + OrigPi0(Length):i 
     + '/	     

c D0
      common /tag3/  
     + myd0Mass(mnum),	
     + myd0kkMCMass(mnum),         myd0piKpMCMass(mnum),
     + myd0piKmMCMass(mnum),       myd0ppMCMass(mnum), 
     + myd0pPpMCMass(mnum),        myd0pPmMCMass(mnum), 
     + myd0kkMass(mnum),         myd0piKpMass(mnum),
     + myd0piKmMass(mnum),       myd0ppMass(mnum), 
     + myd0pPpMass(mnum),        myd0pPmMass(mnum), 
     + myd0kkUPMass(mnum),         myd0piKpUPMass(mnum),
     + myd0piKmUPMass(mnum),       myd0ppUPMass(mnum), 
     + myd0pPpUPMass(mnum),        myd0pPmUPMass(mnum), 
     + myd0VtxFitChisq(mnum),    
     + myd0Pcm(mnum),            myd0FkPromptMass(mnum),
     + myd0fakeMasspKp(mnum),    myd0FakeMasspKm(mnum),
     + myd0swpMass(mnum),        myd0FlightDist(mnum),
     + myd0trueKsPP(mnum),    
     + myd0DecMode(mnum),        myKpPidBit(mnum),
     + myKmPidBit(mnum),         myPipPidBit(mnum),
     + myPimPidBit(mnum),        myTrueD0Flg(mnum),
     + myrecD0CaDec(mnum),       myd0RecDec(mnum)

c D0
      DATA tag3/'  
     + d0Mass(Length):r,	
     + d0kkMCMass(Length):r,         d0piKpMCMass(Length):r,
     + d0piKmMCMass(Length):r,       d0ppMCMass(Length):r, 
     + d0pPpMCMass(Length):r,        d0pPmMCMass(Length):r, 
     + d0kkMass(Length):r,         d0piKpMass(Length):r,
     + d0piKmMass(Length):r,       d0ppMass(Length):r, 
     + d0pPpMass(Length):r,        d0pPmMass(Length):r, 
     + d0kkUPMass(Length):r,         d0piKpUPMass(Length):r,
     + d0piKmUPMass(Length):r,       d0ppUPMass(Length):r, 
     + d0pPpUPMass(Length):r,        d0pPmUPMass(Length):r, 
     + d0VtxFitChisq(Length):r,  
     + d0Pcm(Length):r,            d0FkPromptMass(Length):r,
     + d0fakeMasspKp(Length):r,    d0FakeMasspKm(Length):r,
     + d0swpMass(Length):r,        d0FlightDist(Length):r,
     + d0trueKsPP(Length):i,
     + d0DecMode(Length):i,        KpPidBit(Length):i,
     + KmPidBit(Length):i,         PipPidBit(Length):i,
     + PimPidBit(Length):i,        TrueD0Flg(Length):i,
     + recD0CaDec(Length):i,       d0RecDec(Length):i
     + '/ 

c B
      common /tag4/  
     + myBFCLeo(mnum),           myBFLegendre(mnum), 
     + myCosBThetaT(mnum),       myKbKUpMass(mnum),
     + myKbKlowMass(mnum),       myKbLepMass(mnum),
     + myKbLepMassInt(mnum),     myKbKsUpMass(mnum),      
     + myKbKslowMass(mnum),      myQhemiDiff(mnum),
     + myKroeBvtxDoca(mnum),     mymixNeuMass(mnum),
     + mymixChgMass(mnum),       myKbD0Doca(mnum),
     + myKbD0DocaInLog(mnum),
     + mycosBmomDThrDFm(mnum),   myBCMcosTheta(mnum),  
     + mycosBDVtxDmom(mnum),     myBvtxFitChisq(mnum),
     + myprobB0E(mnum),          myprobB0Mu(mnum), 
     + myprobB0Lep(mnum),        myBandBRoEz(mnum),    
     + myBandBRoeZInLog(mnum),   myBandBRoeY(mnum),
     + mybdVtxVctChisq(mnum),    myexclTruth(mnum),
     + mytrueBflg(mnum),         myDecInChrm(mnum),
     + mymES(mnum),              myDeltaE(mnum),
     + mynnout(mnum),            mybknnout(mnum) 


c B
      DATA tag4/'  
     + BFCLeo(Length):r,           BFLegendre(Length):r, 
     + CosBThetaT(Length):r,       KbKUpMass(Length):r,
     + KbKlowMass(Length):r,
     + KbLepMass(Length):r,        KbLepMassInt(Length):i,
     + KbKsUpMass(Length):r,       
     + KbKslowMass(Length):r,      QhemiDiff(Length):r,
     + KroeBvtxDoca(Length):r,     mixNeuMass(Length):r,
     + mixChgMass(Length):r,       KbD0Doca(Length):r,
     + KbD0DocaInLog(Length):r,
     + cosBmomDThrDFm(Length):r,   BCMcosTheta(Length):r,  
     + cosBDVtxDmom(Length):r,     BvtxFitChisq(Length):r,
     + probB0E(Length):r,          probB0Mu(Length):r, 
     + probB0Lep(Length):i,        BandBRoeZ(Length):r,   
     + BandBRoeZInLog(Length):r,   BandBRoeY(Length):r,
     + bdVtxVctChisq(Length):r,    exclTruth(Length):i,
     + trueBflg(Length):i,         DecInChrm(Length):i,
     + mES(Length):r,              DeltaE(Length):r, 
c neural net work
     + nnout(Length):r,            bknnout(Length):r 
     + '/

      CALL HLIMIT(HSPACE)

      IDOUT=1
c****************D0 K / D0 pi mode*************
       print *, " D0 K / D0 pi mode? "
       print *, " 1 for D0 K, 2 for D0 pi"
       read(*,*) kpimode  

c ********** select the ntuple category *******
       print *, "which sample do you want to select in BB:"
       print *, "0: D0K, 1: D0Pi, 2: other"
       read(*,*) ncat 
c****** dump out the wrong empty ntuple **********
       print *, "dump out the empty ntuple runs?"
       print *,"yes, 1; no, 0"
       read(*,*) dumpfile

c********** out put ntuple *********
       print *, "give a prefix name to the output file,def=file"
       write(*,'(a,$)')
       read(*,'(a)') pref
       if ( pref .eq. ' ' ) pref = 'file'
       call locate(pref,np1,np2)
c**********pi0 mode **********************
c suppress continuum background
           open(unit=84, file='coef_i_3pi.txt',status='old')
            do in3pi=1,13
                read(84,*) 
     >                   hlaycoef(5,in3pi),
     >                   hlaycoef(4,in3pi), hlaycoef(3,in3pi),
     >                   hlaycoef(2,in3pi), hlaycoef(1,in3pi)
           enddo
           close(84)
           open(unit=85,file='coef_o_3pi.txt',status='old')
           do ol3pi=1,13
              read(85,*)  olaycoef(ol3pi)
           enddo
           close(85) 
           open(unit=86,file='act_bias_3pi.txt',status='old')
           do ab=1,19
              read(86,*) act(ab), bias(ab)
           enddo
           close(86)
       print *, " test for coefficient "
       print *, hlaycoef(1,1)," ", hlaycoef(1,2)," ", hlaycoef(1,3)
          

c suppress all background

         open(unit=81, file='coefbk_i_3pi.txt',status='old')
	 tmpnum = 23
         do inkkp=1,tmpnum
                read(81,*) hlaybkcoef(9,inkkp), hlaybkcoef(8,inkkp), 
     >                   hlaybkcoef(7,inkkp), hlaybkcoef(6,inkkp), 
     >		         hlaybkcoef(5,inkkp),
     >                   hlaybkcoef(4,inkkp), hlaybkcoef(3,inkkp),
     >                   hlaybkcoef(2,inkkp), hlaybkcoef(1,inkkp)
           enddo
           close(81)
           open(unit=82,file='coefbk_o_3pi.txt',status='old')
           do olbkkkp=1,tmpnum
              read(82,*)  olaybkcoef(olbkkkp)
           enddo
           close(82)
           open(unit=83,file='act_biasbk_3pi.txt',status='old')
           do ab=1,(10+tmpnum)
              read(83,*) actbk(ab), biasbk(ab)
           enddo
           close(83)
       print *, " test for coefficient "
       print *, hlaycoef(1,1)," ", hlaycoef(1,2)," ", hlaycoef(1,3)


c**********input ntuple************
       if(dumpfile.eq.1) open(unit=89,file='err_run',status='new')
       IDIN=22
       print *, "which file list do u want to input"
       read(5,*) infilelist

      open( unit=90,file=infilelist,status='old')
   
      do 1000 i=1,300000
      	read(90,5001, end=111) infilename
5001    format(A110)

        LREC=0 

       	CALL HROPEN(10,'B2D0KchNonCPAnal',infilename,'P',LREC,ISTAT)
       	if( istat .ne. 0 ) then
       	  write(*,*) 'hropen error ', istat
       	endif
       	CALL HCDIR('//B2D0KchNonCPAnal',' ')
       	ICYCLE=99999
       	IOFF=0

       	CALL HRIN(IDIN,ICYCLE,IOFF)
       	CALL HBNAME(IDIN,' ',0,'$CLEAR')
       	Call HBNAME(IDIN,'EVENT',eventNumber,'$SET')
       	CALL HBNAME(IDIN,'D0KCH',Length,'$SET')
       	CALL HNOENT(IDIN,NCAND)

c******* begin loop over events***********
       	print *, NCAND
        if(NCAND.eq.-1 .and. dumpfile.eq.1) write(89,*) infilename
        if(NCAND .gt. 0 ) print *,infilename
      	do 100 jj=1,NCAND
          CALL HGNT(IDIN,jj,IERR)

	  myLength=0 
          CanIndex = 0

         ncomb=0
         ncomb2=0

         do ii=1, 150
	     index(ii) = 0
         enddo
 
         do  110 ii=1,Length
c pi0 mode 6
           if(d0DecMode(ii).ne.6) goto 110

c CLEO FISHER, ThetaT cos value and D0 K helic
           if(abs(cosBThetaT(ii)).gt.0.8) goto 110
           if(r2.gt.0.5) goto 110
c delta E cut
           if(deltaE(ii).lt.-0.12.or.deltaE(ii).gt.0.14) goto 110  
c HardTrk cuts:
          if(HdTrkThcErr(ii).le.0) goto 110
          if(HdTrkTheta(ii).lt.0.25.or.HdTrkTheta(ii).gt.2.55)
     &          goto 110
          if(HdTrkPcm(ii).lt.0.5) goto 110
          if(HdTrkNPhot(ii).lt.5) goto 110
          if(hdtrkelecpid(ii).gt.0.or.hdtrkmuonpid(ii).gt.0) goto 110 

c d0 daughter trk cuts here:
 	 if(D0DauTrkTheta1(ii).lt.0.25.or.D0DauTrkTheta1(ii).gt.
     &       2.55) goto 110
 	 if(D0DauTrkTheta2(ii).lt.0.25.or.D0DauTrkTheta2(ii).gt.
     &       2.55) goto 110
c         if(D0DauTrkP1(ii).lt.0.25.or.D0DauTrkP2(ii).lt.0.25) goto 110 

c d0 momentum cut here
          if(D0Pcm(ii).lt.1.5) goto 110  

c pid cut
          if(kpimode.eq.1 ) then 
		if(Hdtrkkaonnn(ii)<0.5) goto 110
          endif
          if(kpimode .eq. 2 ) then
		if(Hdtrkkaonnn(ii)>0.5) goto 110
          endif	

c to select the different categories:
c signal B-> D0K 0; charmless 1, d0 pi 2; d0pipi 3; other 4;
        if( ( B1decmode.eq.16.or.B2decmode.eq.16)
     >      .and. decinchrm(ii).ne.0 ) then
                ncatt = 0  ! signal
        else if( ( B1decmode.eq.156.or.B2decmode.eq.156 )
     >           .and. (decinchrm(ii).ne.0)
     >           .and.(B1decmode.ne.16)
     >           .and. (B2decmode.ne.16) ) then
                      ncatt=1  ! dpi
        else
                ncatt = 2
        endif
c        print *, "cat type ", ncatt
        if(ncat .ne. ncatt) goto 110


c pi0 mode only:
   
c B+ -> D0(pi+pi- pi0) K+   
             if(abs(d0mass(ii)-mdpdg).gt.4.5*smd2) goto 110 !d0 mass
             if(abs(pi0mass(ii)-0.135).gt.0.019) goto 110 ! pi0 mass
             if(GamLat1(ii).gt.0.8.or.GamLat2(ii).gt.0.8) goto 110

c            if(pi0cmmom(ii).lt.0.25) goto 110 ! pi0 cm mom
             !!! D0 daughter pid here
             !!! both trks are not kaon 
             if(KpPidBit(ii).ge.19.or.KmPidBit(ii).ge.19) goto 110 
c    how many combinations:
             ncomb2 = ncomb2+1   
             index(ncomb2) = ii    


110     enddo
          
         if(ncomb2.eq.0) goto 100
         rndVal = int(ncomb2*rand(ncomb2)+1.0) 
         CanIndex = index(rndVal)


         ncomb = ncomb2
 
 	
   	   myLength=myLength+1


c neuraon stuff:
           myBandBRoeZInLog(myLength) = log(abs(BandBRoeZ(CanIndex)))  
	   myKbD0docaInLog(myLength) = log(abs(KbD0doca(CanIndex)))

c lepton variables 
           if(abs(probB0Mu(CanIndex)).gt.0.5.and.
     >		abs(probB0e(CanIndex)).gt.0.5) 
     >     	  myprobB0Lep(myLength) = 0 
	   if(abs(probB0Mu(CanIndex)).le.0.1.or.
     >		abs(probB0e(CanIndex)).le.0.1)
     >    	  myprobB0Lep(myLength) = 1
	   if( (abs(probB0Mu(CanIndex)).gt.0.1.and.
     >		abs(probB0Mu(CanIndex)).lt.0.4)
     >    .or.(abs(probB0e(CanIndex)).gt.0.1.and.
     >     	abs(probB0e(CanIndex)).lt.0.4)  )
     >                myprobB0Lep(myLength) = 2
           if( (abs(probB0Mu(CanIndex)).ge.0.4.and.
     > 	       abs(probB0Mu(CanIndex)).le.0.5).or.
     >	       (abs(probB0e(CanIndex)).ge.0.4.and.
     >         abs(probB0e(CanIndex)).le.0.5)) myprobB0Lep(myLength)=3	

	   if(KbLepMass(CanIndex)>1.87) myKbLepMassInt(myLength) = 2
	   if(KbLepMass(CanIndex).le.1.87.and.KbLepMass(CanIndex).gt.0)
     >          myKbLepMassInt(myLength) = 1
           if(KbLepMass(CanIndex).le.0) myKbLepMassInt(myLength) = 0


c dksvtxvctchisq
           if(dksVtxVctChisq(CanIndex).gt.0.)
     >                  mydksVtxVctChisqint(myLength) = 1             
           if(dksVtxVctChisq(CanIndex).gt.5.)
     >                  mydksVtxVctChisqint(myLength) = 2   
           if(dksVtxVctChisq(CanIndex).gt.1000.) 
     >                  mydksVtxVctChisqint(myLength) = 3
           if(dksVtxVctChisq(CanIndex).ge.10000.)
     >                  mydksVtxVctChisqint(myLength) = 4       
           if(dksVtxVctChisq(CanIndex).le.0) 
     >          	mydksVtxVctChisqint(myLength) = 0

c ks decay length
          if(ksDecLen(CanIndex).le.0.)   myKsDecLenInt(myLength) = 0
          if(KsDecLen(CanIndex).gt.0.)   myKsDecLenInt(myLength) = 1
          if(KsDecLen(CanIndex).ge.0.15) myKsDecLenInt(myLength) = 2
          if(KsDecLen(CanIndex).ge.5.)   myKsDecLenInt(myLength) = 3
          if(KsDecLen(CanIndex).ge.20.)  myKsDecLenInt(myLength) = 4

c *******************pi0 mode here**********************************  
c  pi0 mode
c continuum
       do inp = 1, 13
       hlaytmp=0.0
        hlaytmp=BFLegendre(canIndex)*hlaycoef(1,inp) 
     >     +myBandBRoEzInLog(myLength)*hlaycoef(2,inp) 
     >     +CosBThetaT(CanIndex)*hlaycoef(3,inp)
     >     +myKbD0docaInLog(myLength)*hlaycoef(4,inp) 
     >     +myProbB0Lep(myLength)*hlaycoef(5,inp) 

c     >     +myKbLepMassInt(myLength)*hlaycoef(6,inp) 
c     >     +Qhemidiff(CanIndex)*hlaycoef(7,inp) 

        hlaytmp = hlaytmp + bias(5+inp)
        if(hlaytmp.gt.37.) hlaytmp = 37.
        if(hlaytmp.lt.-37) hlaytmp =-37.
          hlay(inp) = 1./(1.0+exp(-hlaytmp))
        hlaytmp = 0 
c      print *, "hid layer node",hlay(inp)
       enddo

        mynnout(myLength) = 0.0
        nntest = 0.0
        do inhl = 1, 13
c     print *,"hid node",inhl," ",hlay(inhl),"  coeff",olaycoef(inhl)
          nntest = nntest + hlay(inhl)*olaycoef(inhl)
        enddo
        nntest = nntest +bias(19)
        if(nntest.gt.37) nntest=37.
        if(nntest.lt.-37) nntest=-37.
        mynnout(myLength) = 1.0/(1.0+exp(-nntest))
c       print *, "nnout ", mynnout(myLength)
        nntest=0


c to suppress all background

              tmpnum=23
              do inp=1,tmpnum
                  hlaytmp = 0.
                  hlaytmp = pi0cmmom(CanIndex)*hlaybkcoef(1,inp)
     >		     + pi0hel(CanIndex)*hlaybkcoef(2,inp)
     >               + pi0mass(CanIndex)*hlaybkcoef(3,inp)
     >               + best2gmamass1(CanIndex)*hlaybkcoef(4,inp)
     >               + best2gmamass2(CanIndex)*hlaybkcoef(5,inp)
     >               + bestcos2gmahel1(CanIndex)*hlaybkcoef(6,inp)
     >               + bestcos2gmahel2(CanIndex)*hlaybkcoef(7,inp)
     >               + cosbmomdthrdfm(CanIndex)*hlaybkcoef(8,inp)
     >               + cosbdvtxdmom(CanIndex)*hlaybkcoef(9,inp)

                   hlaytmp = hlaytmp + biasbk(9+inp)
                   if(hlaytmp.gt. 37.) hlaytmp = 37.
                   if(hlaytmp.lt.-37.) hlaytmp = -37.
                   hlaybk(inp) = 1.0/(1.0+exp(-hlaytmp))
                   hlaytmp=0.

               enddo
               nntest=0
               do inp=1,tmpnum
                  nntest= nntest+hlaybk(inp)*olaybkcoef(inp)
               enddo
               nntest = nntest + biasbk(tmpnum+10)
               if(nntest.gt. 37.) nntest= 37.
               if(nntest.lt.-37.) nntest=-37.
               mybknnout(myLength) = 0.0
               mybknnout(myLength) = 1.0/(1.0+exp(-nntest))
               nntest = 0

c===============================================================
c kaon  	   
           myHdTrkThc(myLength)= HdTrkThc(CanIndex) 
           myHdTrkThcErr(myLength)= HdTrkThcErr(CanIndex) 
           myHdTrkKaonNN(mylength)= HdTrkKaonNN(CanIndex)
           myHdTrkCMTheta(mylength)= HdTrkCMTheta(CanIndex)
           myHdTrkPcm(mylength)= HdTrkPcm(CanIndex)
           myHdTrkChge(mylength)= HdTrkChge(CanIndex)
           myHdtrkTruth(mylength) = HdtrkTruth(CanIndex)
c Ks
           myKsMass(mylength)= KsMass(CanIndex)
           myKsDecLen(mylength)= KsDecLen(CanIndex)
           myKsVtxFitChisq(mylength) = KsVtxFitChisq(CanIndex)
           myKsHel(mylength) = KsHel(CanIndex)
           mycosDKsVtx(mylength)= cosDKsVtx(CanIndex)
           mydksVtxVctChisq(mylength)= dksVtxVctChisq(CanIndex)
           myOrigKs(mylength) = OrigKs(CanIndex)
c pi0
           mypi0Mass(mylength)= pi0Mass(CanIndex)
           mypi0Hel(mylength) = pi0Hel(CanIndex)
           mypi0CMmom(mylength)= pi0CMmom(CanIndex)
           mypi0EAsy(mylength) = pi0EAsy(CanIndex)
	   myBest2GmaMass1(mylength)=Best2GmaMass1(CanIndex)
	   myBest2GmaMass2(mylength)=Best2GmaMass2(CanIndex)
	   myBestcos2gmaHel1(mylength)=Bestcos2gmaHel1(CanIndex)
	   myBestcos2gmaHel2(mylength)=Bestcos2gmaHel2(CanIndex)
	   myBest2gPcm1(mylength) = Best2gPcm1(CanIndex)
	   myBest2gPcm2(mylength) = Best2gPcm2(CanIndex)
 	   myOrigPi0(mylength) = OrigPi0(CanIndex)
c D0
           myd0Mass(mylength)= d0Mass(CanIndex)
           myd0kkMCMass(mylength)= d0kkMCMass(CanIndex)
           myd0piKpMCMass(mylength)= d0piKpMCMass(CanIndex)
           myd0piKmMCMass(mylength)= d0piKmMCMass(CanIndex)
           myd0ppMCMass(mylength)= d0ppMCMass(CanIndex)
           myd0pPpMCMass(mylength)= d0pPpMCMass(CanIndex)
           myd0pPmMCMass(mylength)= d0pPmMCMass(CanIndex)
           myd0kkMass(mylength)= d0kkMass(CanIndex)
           myd0piKpMass(mylength)= d0piKpMass(CanIndex)
           myd0piKmMass(mylength)= d0piKmMass(CanIndex)
           myd0ppMass(mylength)= d0ppMass(CanIndex)
           myd0pPpMass(mylength)= d0pPpMass(CanIndex)
           myd0pPmMass(mylength)= d0pPmMass(CanIndex)
           myd0kkUPMass(mylength)= d0kkUPMass(CanIndex)
           myd0piKpUPMass(mylength)= d0piKpUPMass(CanIndex)
           myd0piKmUPMass(mylength)= d0piKmUPMass(CanIndex)
           myd0ppUPMass(mylength)= d0ppUPMass(CanIndex)
           myd0pPpUPMass(mylength)= d0pPpUPMass(CanIndex)
           myd0pPmUPMass(mylength)= d0pPmUPMass(CanIndex)
           myd0VtxFitChisq(mylength)= d0VtxFitChisq(CanIndex)
           myd0Pcm(mylength)= d0Pcm(CanIndex)
           myTrueD0Flg(mylength) = TrueD0Flg(CanIndex) 
           myd0FkPromptMass(mylength) = d0FkPromptMass(CanIndex)
           myd0fakeMasspKp(mylength)=d0fakeMasspKp(CanIndex)
           myd0fakeMasspKm(mylength)=d0fakeMasspKm(CanIndex)
           myd0swpMass(mylength) = d0swpMass(CanIndex)
	   myd0flightDist(mylength) = d0flightDist(CanIndex)
           myd0trueKsPP(mylength) = d0trueKsPP(CanIndex)
           myd0DecMode(mylength)= d0DecMode(CanIndex)
           myKpPidBit(mylength)= KpPidBit(CanIndex)
           myKmPidBit(mylength)= KmPidBit(CanIndex)
           myPipPidBit(mylength)= PipPidBit(CanIndex)
           myPimPidBit(mylength)= PimPidBit(CanIndex)
           mybdVtxVctChisq(mylength)= bdVtxVctChisq(CanIndex)
	   myrecD0CaDec(mylength) = recD0CaDec(CanIndex)
	   myd0RecDec(mylength) = d0RecDec(CanIndex)		
c B
           myBFCLEO(mylength)= BFCLEO(CanIndex)
           myBFLegendre(mylength)= BFLegendre(CanIndex)
           myCosBThetaT(mylength)= CosBThetaT(CanIndex)
	   myKbKUpMass(mylength)= KbKUpMass(CanIndex)
	   myKbKlowMass(mylength)= KbKlowMass(CanIndex)	
	   myKbKsUpMass(mylength)= KbKsUpMass(CanIndex)
	   myKbKslowMass(mylength)= KbKslowMass(CanIndex)
	   myKbLepMass(mylength)= KbLepMass(CanIndex)
           myQhemiDiff(mylength)= QhemiDiff(CanIndex)
           myKroeBvtxDoca(mylength)= kroeBvtxDoca(CanIndex)
	   mymixNeuMass(mylength)= mixNeuMass(CanIndex)
	   myMixChgMass(mylength)= MixChgMass(CanIndex)
           myKbD0Doca(mylength)= KbD0Doca(CanIndex)
	   mycosBmomDThrDFm(mylength)= cosBmomDThrDFm(CanIndex)
           myBCMcosTheta(mylength)= BCMcosTheta(CanIndex)
           mycosBDVtxDmom(mylength)= cosBDVtxDMom(CanIndex)
           myBvtxFitChisq(mylength)= BvtxFitChisq(CanIndex)
	   myprobB0E(mylength)= probB0E(CanIndex)
	   myProbB0Mu(mylength)= ProbB0Mu(CanIndex)
           myBandBRoeZ(mylength)= BandBRoeZ(CanIndex)
           myBandBRoeY(mylength)= BandBRoeY(CanIndex)
           mybdVtxVctChisq(mylength)= bdVtxVctChisq(CanIndex)	
           myexclTruth(mylength)= exclTruth(CanIndex)
           mytrueBflg(mylength)= trueBflg(CanIndex)
           myDecInchrm(mylength)= DecInchrm(CanIndex)
           mymES(mylength)= mES(CanIndex)
           myDeltaE(mylength)= DeltaE(CanIndex)


         myeventNumber=eventNumber
	 myrunNumber= runNumber 
	 myLower= Lower
 	 myUpper= Upper
	 myR2= r2
	 myB1decmode = B1decmode
	 myB2decMode = B2decMode
         myD1CaDecMode = D1CaDecMode
         myD2CaDecMode = D2CaDecMode
	 myHem1Mass = Hem1Mass
	 myHem2Mass = Hem2Mass	

***********************
*start to write output
***********************
      nwrite=mod(iwrt,15000)
      nfile=1+iwrt/15000

      if(nwrite.eq.0.and.iwrt.ne.0) then
     	Call hcdir('//NTNEW',' ')
      	call HROUT(0,ICYCLE,' ')
      	call HREND('NTNEW')
      	close (31)
      endif

      iwrt=iwrt+1

      if(nwrite.eq.0) then
      	write(fname,*) nfile,'.hbk'
      	print *,'output filename = ',pref(np1:np2)//"2-"//fname(2:)
     	LREC=8190
      	ISTAT = 0
       	call hropen(31,'NTNEW',pref(np1:np2)//"2-"//fname(2:),
     &     'N', LREC, ISTAT)
       	call hcdir('//NTNEW',' ')
       	call hbnt(idout,'datantple',' ')
       	call hbname(idout,'Event',myeventNumber,tag1)
       	call hbname(idout,'D0Kch',myLength,tag2)
       	call hbname(idout,'D0Kch',myD0Mass(1),tag3)
       	call hbname(idout,'D0Kch',myBFCLeo(1),tag4)
      endif

c************fill in the new ntuple********
        call hcdir('//NTNEW',' ')
        call hfnt(idout)

100    continue

c******* close current input ntuple*******
      CALL HCDIR('//B2D0KchNonCPAnal',' ')
      call hrend('B2D0KchNonCPAnal')
      close(10)

1000  continue
111   continue
      close(90)

      Call hcdir('//NTNEW',' ')
      call HROUT(0,ICYCLE,' ')
      call HREND('NTNEW')
      close (31)
      if(dumpfile.eq.1) close(89)
      END

        subroutine locate(chars,n1,n2)
        IMPLICIT NONE
        character*255 chars
        integer n1, n2
        if ( chars.eq.' ' ) then
          n1=0
          n2=0
          return
        endif
        do n1=1,255
          if ( chars(n1:n1).ne.' ' ) goto 100
        enddo
100     do n2=n1+1,255
          if ( chars(n2:n2).eq.' ' ) goto 200
        enddo
200     n2=n2-1
        end

