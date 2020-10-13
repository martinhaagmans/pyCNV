import os

HOME = '/data/dnadiag'

dbdir = os.path.join(HOME,  'databases', 'cnv')
outputdir = os.path.join(HOME,  'cnvoutput')
pipelinedir = os.path.join(HOME, 'ngstargets')

patientinfo = os.path.join(HOME,  'databases', 'patientinfo.sqlite')

dbgeneral = 'general.sqlite'
poscontable = 'poscontrols'
badsampletable = 'badsamples'
badregiontable = 'toexclude'

