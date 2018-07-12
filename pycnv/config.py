import os

HOME = os.path.expanduser('~')

dbdir = os.path.join(HOME, 'Documents', 'cnv', 'databases')
outputdir = os.path.join(HOME, 'Documents', 'cnv', 'output')
pipelinedir = os.path.join(HOME, 'Documents', 'ngstargets')

dbgeneral = 'general.sqlite'
poscontable = 'poscontrols'
badsampletable = 'badsamples'
badregiontable = 'toexclude'
