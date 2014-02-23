# source this file in (t)csh to create all three Dalitz projections

nice +5 bbrroot -b -q setup.cc 'plotDalitz.cc(m12)' > & ! plotDalitz-m12.log &
nice +5 bbrroot -b -q setup.cc 'plotDalitz.cc(m13)' > & ! plotDalitz-m13.log &
nice +5 bbrroot -b -q setup.cc 'plotDalitz.cc(m23)' > & ! plotDalitz-m23.log &
