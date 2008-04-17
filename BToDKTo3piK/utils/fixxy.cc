// $Id: fixxy.cc,v 1.1 2006/05/08 23:13:19 fwinkl Exp $
// Fix x and y

void fixxy(Bool_t fix = kTRUE)
{
  dalitzHolderP.sigGoodD0Type().x()->setConstant(fix);
  dalitzHolderP.sigGoodD0Type().y()->setConstant(fix);
  dalitzHolderN.sigGoodD0Type().x()->setConstant(fix);
  dalitzHolderN.sigGoodD0Type().y()->setConstant(fix);
}
