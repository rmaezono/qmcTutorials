#!/usr/bin/python

from sys import exit

interactive = False

class gpifParamError(Exception):
  def __init__(self,param):
    self.param = param

def usage():
### Prints a message showing the correct usage of the script
  print '''
Correct usage: plot.py [OPTIONS] datafiles

OPTIONS from the following:

-h, --help    : display this message

--interactive : open with interactive mode. Allows interactive correction of
                errors in input file

    -i <file> : specify an input file. Default is 'input'. If the file does not
                exist, an example is copied from LOUIS examples directory

    -o <file> : output for plots to file. Ignored for x11 or wxt output.
                File extension is decided automatically

    -l <file> : read location of data files from file. Expects one location per
                line

    -s <file> : save the gnuplot batch file to the location specified.

    -r [file] : replots from batch file. If no file provided, uses '.gpout.gp'
                which will always be the last plot made in the directory.
'''


def cont():
### Get yes/no response from user. Returns True or False accordingly
  cont = ''
  while cont not in ('y', 'yes', 'n', 'no'):
    cont = raw_input('>>> ')
    if cont in ('y', 'yes'):
      return True
    elif cont in ('n', 'no'):
      return False


def datatype(datafile):
### Attempt to recover type of data in file by assuming a sensible naming
### convention....
  if datafile.find('trajectory') != -1:
    return 'trajectory'
  elif datafile.find('raw') != -1:
    return 'raw'
  elif datafile.find('cg') != -1:
    return 'coarse-grained'
  elif datafile.find('smoothed') != -1:
    return 'smoothed'
  elif datafile.find('h') != -1:
    return 'coarse-grained'
  else:
    return 'unknown'


def setview(gpif):
### Get appearance options from input file and set in gnuplot
  try:
    colour = gpif.get('colour','boolean')
  except gpifParamError as err1:
    try:
      colour = gpif.get('color','boolean')
    except gpifParamError as err2:
      print 'Neither \'colour\' nor \'color\' found in input file'
      exit()
  autorotate = gpif.get('autorotate', 'boolean')
  try:
    key = gpif.get('colour_key','boolean')
  except gpifParamError as err1:
    try:
      key = gpif.get('color_key','boolean')
    except gpifParamError as err2:
      print 'Neither \'colour_key\' nor \'color_key\' found in input file'
      exit()
  dim = gpif.get('dimensionality','integer')
  plot_type = gpif.get('plot_type')
  viewstr = ''
  if dim == 2 and plot_type in ('map','trajectory'):
    viewstr += 'set size ratio -1\n'
  if not colour:
    viewstr += 'set palette defined (0 0.1 0.1 0.1, 0.5 0.7 0.7 0.7,\
 1 0.9 0.9 0.9)\n'
  viewstr += 'set xyplane 0\n'
  if not autorotate:
    rot_x = gpif.get('rot_x')
    rot_z = gpif.get('rot_z')
    viewstr += 'set view '+rot_x+','+rot_z+'\n'
  if not key:
    viewstr += 'unset colorbox\n'
  return viewstr


def setterm(term,output,colour,n_plot=0,fps=0):
### Set the terminal according to plot_ouput in the input file
  if colour:
    colourlist = 'xffffff x000000 x000000 xc00000 x00c000 x0000c0 xc0c000\
 xc000c0 x00c0c0 x000000'
  else:
    colourlist = 'xffffff x000000 x000000 x202020 x404040 x606060 x808080\
 xa0a0a0 xc0c0c0 x000000'

  termstr = ''

  if term in ('wxt', 'x11'):
    termstr += 'set terminal '+term+' \''+str(n_plot)+'\' enhanced\n'

  elif term in ('postscript','ps'):
    if n_plot == 0:
      termstr += 'set terminal postscript enhanced '+('colour solid' if colour
        else 'monochrome')+' linewidth 0.5\n'
    termstr += 'set output \''+output+'.ps\'\n'
    print 'Plot will be saved as \''+output+'.ps\''

  elif term == 'eps':
    if n_plot == 0:
      termstr += 'set terminal postscript eps enhanced '+('colour solid' if
        colour else 'monochrome')+' linewidth 0.5\n'
    termstr += 'set output \''+output+'.eps\'\n'
    print 'Plot will be saved as \''+output+'.eps\''

  elif term in ('png','gif'):
    if n_plot == 0:
      termstr += 'set terminal '+term+' size 720,720 enhanced '+colourlist+'\n'
    termstr += 'set output \''+output+'.'+term+'\n'
    print 'Plot will be saved as \''+output+'.'+term+'\''
    return termstr

  elif term in ('jpg','jpeg'):
    if n_plot == 0:
      termstr += 'set terminal jpeg size 720,720 enhanced '+colourlist+'\n'
    termstr += 'set output \''+output+'.jpg\n'
    print 'Plot will be saved as \''+output+'.jpg\''

  elif term == 'agif':
    delay_time = 100.0/fps
    termstr += 'set terminal gif transparent animate delay '+str(delay_time)+\
' optimize size 720,720 enhanced '+colourlist+'\n'
    termstr += 'set output \''+output+'.gif\'\n'
    print 'Plot will be saved as \''+output+'.gif\''

  elif term == 'svg':
    if n_plot == 0:
      termstr += 'set terminal svg enhanced\n'
    termstr += 'set output \''+output+'.svg\'\n'
    print 'Plot will be saved as \''+output+'.svg\''

  elif term in ('pdf','pdfcairo'):
    if n_plot == 0:
      termstr += 'set terminal pdfcairo enhanced '+('colour solid' if colour
        else 'monochrome')+' linewidth 0.5\n'
    termstr += 'set output \''+output+'.pdf\'\n'
    print 'Plot will be saved as \''+output+'.pdf\''
  
  else:
    termstr += 'set terminal '+term+'\n'
    print 'Output to file? (y\\n)'
    if cont():
      print 'Enter name of output file'
      output = ''
      while output == '':
        output = raw_input('>>> ')
      termstr += 'set output \''+output+'\'\n'

  return termstr


def settitle(opt,datafile):
  if opt == 'none':
    return 'unset title\n'

  elif opt == 'prompt':
    print '\nPlease enter plot title for file \''+datafile+'\':'
    titlestring = raw_input('>>> ')
    return 'set title \''+titlestring+'\'\n'

  elif opt == 'auto':            # Return a sensible title for the plot if name of datafile is in default format for LOUIS output. Otherwise unset title
    if 'density' in datafile:
      title = 'particle density'
      if 'smoothed' in datafile:
        title = 'smoothed '+title
      elif 'cg' in datafile:
        title = 'coarse-grained '+title
      elif 'raw' in datafile:
        title = 'raw '+title
      parts = datafile.split('_')
      for part in parts:
        if 't=' in part:
          time = part.lstrip('t=').rstrip('.dat')
          title += ' at t = '+time
        if 'cg=' in part:
          cglength = part.lstrip('cg=').rstrip('.dat')
          title += ': coarse-graining length '+cglength
      titlestr = 'set title \''+title+'\'\n'
    elif 'psisq' in datafile:
      title = '|psi|^2'
      if 'smoothed' in datafile:
        title = 'smoothed '+title
      elif 'cg' in datafile:
        title = 'coarse-grained '+title
      elif 'raw' in datafile:
        title = 'raw '+title
      parts = datafile.split('_')
      for part in parts:
        if 't=' in part:
          time = part.lstrip('t=').rstrip('.dat')
          title += ' at t = '+time
        if 'cg=' in part:
          cglength = part.lstrip('cg=').rstrip('.dat')
          title += ': coarse-graining length '+cglength
      titlestr = 'set title \''+title+'\'\n'
    elif 'h_integrand' in datafile:
      title = 'integrand of H-function'
      parts = datafile.split('_')
      for part in parts:
        if 't=' in part:
          time = part.lstrip('t=').rstrip('.dat')
          title += ' at t = '+time
        if 'cg=' in part:
          cglength = part.lstrip('cg=').rstrip('.dat')
          title += ': coarse-graining length '+cglength
      titlestr = 'set title \''+title+'\'\n'
    elif 'trajectory' in datafile:
      title = 'trajectory of particle'
      titlestr = 'set title \''+title+'\'\n'
    else:
      print 'Automatic title failed'
      titlestr = 'unset title'
    return titlestr

  else:
    print 'Title option not recognised: interpreting as a string'
    return 'set title \''+opt+'\'\n'







