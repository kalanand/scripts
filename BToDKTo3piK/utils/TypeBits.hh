// $Id: TypeBits.hh,v 1.2 2006/04/08 04:24:24 fwinkl Exp $
// Definition of event type bits

enum {
  SIG_G_BIT = 1,
  SIG_B_BIT = 2,
  DPI_G_BIT = 4,
  DPI_B_BIT = 8,
  DPIX_BIT  = 16,
  DKX_BIT   = 32,
  BB_G_BIT  = 64,
  BB_B_BIT  = 128,
  QQ_G_BIT  = 256,
  QQ_B_BIT  = 1024,

  ALL_TYPES = SIG_G_BIT+SIG_B_BIT+DPI_G_BIT+DPI_B_BIT+DPIX_BIT+DKX_BIT+
              BB_G_BIT+BB_B_BIT+QQ_G_BIT+QQ_B_BIT
};

TString printBits(int i){
  TString result;
  if (i & SIG_G_BIT) result += "SIG_G_BIT ";
  if (i & SIG_B_BIT) result += "SIG_B_BIT ";
  if (i & DPI_G_BIT) result += "DPI_G_BIT ";
  if (i & DPI_B_BIT) result += "DPI_B_BIT ";
  if (i & DPIX_BIT)  result += "DPIX_BIT ";
  if (i & DKX_BIT)   result += "DKX_BIT ";
  if (i & BB_G_BIT) result += "BB_G_BIT ";
  if (i & BB_B_BIT) result += "BB_B_BIT ";
  if (i & QQ_G_BIT) result += "QQ_G_BIT ";
  if (i & QQ_B_BIT) result += "QQ_B_BIT ";
  return result;
}

