/*****************************************************************************
 * Project: RooFit                                                           
 * Package: BToDKTo3piK                                                     
 *    File: $Id: BdkSplitData.cc,v 1.1 2006/04/14 21:17:16 fwinkl Exp $
 * Authors:                                                                  
 *   Frank Winklmeier, Colorado State University                             
 * Description:                                                              
 *   Class to split a RooAbsData into several chunks.
 *                                                                           
 * Copyright (c) 2006, Colorado State University                             
 *****************************************************************************/

#include "BToDKTo3piK/BdkSplitData.hh"

#include "RooFitCore/RooAbsData.hh"
#include "RooFitCore/RooRandom.hh"

ClassImp(BdkSplitData)


/// Constructor
BdkSplitData::BdkSplitData(const char *name, const char *title,
                           RooAbsData& data) :
  TNamed(name, title),
  _chunkSize(0),
  _nChunks(0),
  _thisChunk(0),
  _data(&data)
{
}


/// Copy constructor
BdkSplitData::BdkSplitData(const BdkSplitData& other, const char* name) :
  TNamed(other),
  _chunkSize(other._chunkSize),
  _nChunks(other._nChunks),
  _thisChunk(other._thisChunk),
  _data(other._data)
{
  if (name) SetName(name);
}


/// Set the number of chunks by chunk size
Bool_t BdkSplitData::setChunkSize(Int_t entries)
{
  if (entries <= 0) return kFALSE;

  Int_t N = _data->numEntries();
  if (entries > N) return kFALSE;

  _chunkSize = entries;
  _nChunks = N/_chunkSize;
  resetChunks();

  return kTRUE;
}

/// Set the number of chunks
Bool_t BdkSplitData::setNumChunks(Int_t numChunks)
{
  if (numChunks<=0) return kFALSE;

  Int_t size = _data->numEntries()/numChunks;
  return setChunkSize(size);
}

/// Return chunk number n
/// Caller owns dataset
RooAbsData* BdkSplitData::chunk(Int_t n)
{
  if (n<0 || n>=_nChunks) return 0;
  if (_chunkSize<=0) return 0;

  Int_t start = _chunkSize*n;

  RooAbsData* chunk = _data->emptyClone();
  for (Int_t i = start; i<start+_chunkSize; i++) {
    chunk->add(*_data->get(i));
  }
  
  return chunk;
}
 
/// Return the next chunk (0 if none available)
/// Caller owns dataset
RooAbsData* BdkSplitData::nextChunk()
{
  if (_thisChunk < _nChunks) return chunk(_thisChunk++);
  else return 0;
}


/// Get a random chunk
/// Caller owns dataset
RooAbsData* BdkSplitData::randomChunk()
{
  Int_t n = RooRandom::integer(_nChunks);
  return chunk(n);  
}


/// Rewind to the first chunk
void BdkSplitData::resetChunks()
{
  _thisChunk = 0;
}
