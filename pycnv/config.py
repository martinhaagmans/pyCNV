import os

HOME = os.path.expanduser('~')

dbdir = os.path.join(HOME, 'Documents', 'databases', 'cnv')
outputdir = os.path.join(HOME, 'Documents', 'cnvoutput')
pipelinedir = os.path.join(HOME, 'Documents', 'ngstargets')

dbgeneral = 'general.sqlite'
poscontable = 'poscontrols'
badsampletable = 'badsamples'
badregiontable = 'toexclude'
