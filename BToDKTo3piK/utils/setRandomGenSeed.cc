UInt_t setRandomGenSeed(const UInt_t randomSeed) 
{
  RooRandom::randomGenerator()->SetSeed(randomSeed);
  return randomSeed;
}
