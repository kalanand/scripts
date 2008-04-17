sourceFoundFile GeneratorsQA/EvtGen.tcl


mod talk GfiEvtGen
  GENERATE set continuum
  transformEvent set false

  UDECAY set ../BToDKTo3piK/generator/257_ABAC190_ACAB0_BCAC-4_DDbar.dec
#  UDECAY set ../BToDKTo3piK/generator/257_ABAC190_ACAB0_BCAC-4.dec
#  UDECAY set ../BToDKTo3piK/generator/257_ABAC170_ACAB0_BCAC-4.dec
#  UDECAY set ../BToDKTo3piK/generator/77_ABAC10_ACAB0_BCAC-4.dec
#  UDECAY set ../BToDKTo3piK/generator/77_ABAC10_ACAB0_CBCA-4.dec
#  UDECAY set ../BToDKTo3piK/generator/77_ABAC170_ACAB0_BCAC-4.dec
exit

module talk HbkTupleEnv
  histFileName set 257_ABAC190_ACAB0_BCAC-4_DDbar.hbk
#  histFileName set 257_ABAC190_ACAB0_BCAC-4.hbk
#  histFileName set 257_ABAC170_ACAB0_BCAC-4.hbk
#  histFileName set 77_ABAC10_ACAB0_BCAC-4.hbk
#  histFileName set 77_ABAC10_ACAB0_CBCA-4.hbk
#  histFileName set 77_ABAC170_ACAB0_BCAC-4.hbk
exit

mod talk GefSelectFilter
  BooNew D = GefPdtList
  BooObjects D or D0

  BooNew 2D = GefPdtList
  BooObjects  2D and D0 D0
  BooNew no2D = !2D

  BooCompose Dfilt = and no2D D


  BooNew anti-D = GefPdtList
  BooObjects anti-D or anti-D0

  BooNew 2anti-D = GefPdtList
  BooObjects  2anti-D and anti-D0 anti-D0
  BooNew no2anti-D = !2anti-D

  BooCompose anti-Dfilt = and no2anti-D anti-D


  BooCompose filt = or anti-Dfilt Dfilt

  beforeFilter set filt
exit


mod talk RandomControl
maxEventsPerRun set 50000
exit

#mod talk BtaLoadMcCandidates
#  numEventsToPrint set 0
#exit

mod talk D0
  printTerminals set ""
exit

mod talk anti_D0
  printTerminals set ""
exit

ev begin -nev 50000
exit

