import FWCore.ParameterSet.Config as cms
process = cms.Process("PICKEVENTS")

process.source = cms.Source ("PoolSource",
                                     fileNames = cms.untracked.vstring(
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/FCEB503A-7491-DF11-961F-003048F024F6.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/F69A46A0-8391-DF11-950F-003048F118E0.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/ECAECDE4-7B91-DF11-B97F-001617C3B77C.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/E2027DB3-7791-DF11-B264-003048F1182E.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/D223F468-7191-DF11-B564-003048F118C4.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/CC5A9831-7491-DF11-BC68-001617DBCF6A.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/C61D6098-5B91-DF11-809B-003048CFB40C.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/AA06A857-7D91-DF11-A701-0019B9F72CE5.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/9E8C6DFF-6391-DF11-A73B-001617C3B6DE.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/9A094321-8091-DF11-836D-003048F1C58C.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/8C98E8AC-6991-DF11-9B73-001D09F29849.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/8A790879-6C91-DF11-BACF-001D09F28EC1.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/78A4DC68-7191-DF11-866B-003048F11C5C.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/70433E18-7991-DF11-B9EA-0030487A18A4.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/6CDE8454-7D91-DF11-BE29-001D09F25456.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/4CF657D3-7991-DF11-8D7A-001D09F2438A.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/4C20F2CE-7391-DF11-8447-0030487CD76A.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/4C101C68-7F91-DF11-87C7-000423D33970.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/48DD1D95-6E91-DF11-9E8D-001D09F2447F.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/36A96B1F-8791-DF11-A3A1-003048F0258C.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/30539D39-7491-DF11-B5C1-001617C3B65A.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/1A66FBA6-8391-DF11-AE79-003048F118DE.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/1831AF9B-6C91-DF11-8A0E-0019B9F709A4.root',
    '/store/data/Run2010A/EG/RECO/v4/000/140/331/163F1611-6B91-DF11-A342-003048D2C092.root'
    ),
    eventsToProcess = cms.untracked.VEventRange('140331:384695957')
  )

process.Out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string('HigmMass_ee.root')
                               )

process.e = cms.EndPath(process.Out)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
