import sys
from plot_utils import gpifParamError as ParamError

class inputFile:           # Class to parse input file and return requested
                           # parameters
 def __init__(self, path, interactive=False):
  self.interactive = interactive
  fInput = open(path,'r')  # Open input file
  strInput = []
  self.optarg = []
  for line in fInput:
   if line.strip() == '':    # Blank line
    pass
   elif line.strip()[0] == '#':      # Comment line
    pass
   elif line.strip()[0] == '%':      # Start of an esdf block
    line = fInput.next().strip()
    while line[0] != '%':
     line = fInput.next().strip()
   else:                     # Remove inline comments and add to list of commands
    strInput.append(str(line).split('#')[0].lower())
  for line in strInput:
   if line[0] in ('%',' '):  # Ignore block data in section parsed by esdf
    del line
   else:                   # Add to list of (opt,arg) tuples
    try:
     self.optarg.append((line.split(':')[0].strip(),line.split(':')[1].strip()))
    except IndexError:
     print 'Warning: gpif unable to parse line:\n'+line

 def get(self, param, rtnType='string'):
  if rtnType in ('string', 'str', 's'):
   rtnType = 'str'
   for opt, arg in self.optarg:
    if opt == param:
     return arg
#   return self.paramError(param)
   raise ParamError(param)
  elif rtnType in ('bool', 'boolean', 'b'):
   rtnType = 'bool'
   for opt, arg in self.optarg:
    if opt == param:
     if arg in ('t', 'true', '.true.', '1'):
      return True
     elif arg in ('f', 'false', '.false.', '0'):
      return False
     else:
      return self.valueError(arg, param, rtnType)
#   return self.paramError(param, 'boolean') 
   raise ParamError(param)
  elif rtnType in ('int', 'integer', 'i'):
   rtnType = 'int'
   for opt, arg in self.optarg:
    if opt == param:
     if arg.isdigit():
      return int(arg)
     else:
      return self.valueError(arg, param, rtnType)
#   return self.paramError(param, 'integer')
   raise ParamError(param)
  elif rtnType in ('float', 'f'):
   rtnType = 'float'
   for opt, arg in self.optarg:
    if opt == param:
     try:
      return float(arg)
     except ValueError:
      return self.valueError(arg, param, rtnType)
#   return self.paramError(param, 'float')      
   raise ParamError(param)
  else:
   self.typeError(param, rtnType) 

 def paramError(self, param, rtnType='string'):
  print 'Parameter \''+param+'\' not found in input file.'
  opt = None
  if self.interactive:
   while not opt in ('y', 'n', 'yes', 'no'):
    opt = raw_input('Enter value manually? (y\\n): ')
    if opt in ('y', 'yes'):
     if rtnType == 'string':
      return str(raw_input('>>> '))
     elif rtnType == 'boolean':
      print 'Note: boolean type required for '+param
      rtn = None
      while rtn not in ('t', 'true', 'f', 'false'):
       rtn = raw_input('>>> ').lower()
       if rtn in ('t', 'true'):
        return True
       elif rtn in ('f', 'false'):
        return False
     elif rtnType == 'integer':
      print'Note: integer type required for '+param
      rtn = ''
      while not rtn.isdigit():
       rtn = raw_input('>>> ')
      return int(rtn)
     elif rtnType == 'float':
      print 'Note: float type required for '+param
      rtn = None
      while rtn == None:
       try:
        rtn = float(raw_input('>>> '))
       except ValueError:
        rtn = None
     return rtn
    elif opt in ('n', 'no'):
     sys.exit()
  else:
   sys.exit()

 def typeError(self, param, rtnType):
  print 'Type \''+rtnType+'\' specified for parameter \''+param+'\' not recognised.'
  print 'Currently supported types:\nstring\nboolean\n'
  sys.exit()

 def valueError(self, arg, param, rtnType):
  print '\''+arg+'\' not valid argument of type '+rtnType+' for parameter \''+param+'\''
  if self.interactive:
   opt = None
   while not opt in ('y', 'yes', 'n', 'no'):
    opt = raw_input('Enter a value manually? (y\\n): ')
    if opt in ('y', 'yes'):
     if rtnType == 'bool':
      print 'Please enter an argument of boolean type'
      while not arg in ('t', 'true', 'f', 'false'):
       arg = raw_input('>>> ')
       if arg in ('t', 'true'):
        return True
       elif arg in ('f', 'false'):
        return False
     elif rtnType == 'str':
      print 'Please enter an argument of integer type'
      while not arg.isdigit():
       arg = raw_input('>>> ')
      return int(arg)
     else:
      print 'Unsupported arument type'
      sys.exit()
    elif opt in ('n', 'no'):
     sys.exit()
  else:
   sys.exit()



