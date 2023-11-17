#!/usr/bin/python
#------------------PLOT_SINGLE.PY------------------#
# Define functions that return the string of       #
# gnuplot commands to plot the different types of  #
# LOUIS data files                                 #
#--------------------------------------------------#

import os.path

# ----------------------------- 1D PLOT OPTIONS ----------------------------- #

def plot_1d_trajectory(datafile,n_plot=0):
# Look through file to find number of trajectories
  raw_data = open(datafile, 'r')
  traj_no = 1
  for line in raw_data:
    if line.strip() == '':     # blank lines divide different trajectories
      traj_no += 1
  raw_data.close()
  print str(traj_no)+' trajectories detected'

  plotstr = 'plot \''+datafile+'\' every :::0::0 with line ls 1 notitle'
  if traj_no > 1:              # more than one trajectory
    for i in range(1,traj_no):
      plotstr += ',\\\n \''+datafile+'\' every :::'+str(i)+'::'+str(i)+' with\
 line ls '+str(i)+' notitle'
  plotstr += '\n'
  return plotstr


def plot_1d_density(datafile,n_plot=0):
  plotstr = 'plot '+'\''+datafile+'\' with lines ls 1 notitle\n'
  return plotstr


def plot_1d_bar(datafile,n_plot=0):
  try:
    bar_datafile = datafile.split('_')[0]+'_bar_'+datafile.split('_')[2]
  except IndexError:
    bar_datafile = 'bar_data_'+str(n_plot)+'.dat'
# Generate file containing bar data
##### DO THIS 'IN PLACE' AS IN THE 2D CASE
  raw_data = []
  raw_datafile = open(datafile, 'r')
  for line in raw_datafile:
    if line.strip() == '':   # ignore any blank lines
      pass
    elif line[0] == '#':     # ignore any comment lines
      pass
    else:
      raw_data.append([float(x) for x in line.strip().split()])
  raw_datafile.close()
  delta_x = 0.5 * (raw_data[1][0] - raw_data[0][0])
  plot_data = []
  for point in raw_data:
    plot_data.append([point[0]-delta_x,point[1]])
    plot_data.append([point[0]+delta_x,point[1]])
# Output to file
  plotfile = open(bar_datafile, 'w')
  for point in plot_data:
    plotfile.write(str(point[0])+' '+str(point[1])+'\n')
  plotfile.close()
# Create bar\step type plot
  plotstr = 'plot \''+bar_datafile+'\' with lines ls 1 notitle\n'
  return plotstr

# ----------------------------- 2D PLOT OPTIONS ----------------------------- #

def plot_2d_trajectory(datafile,spacetime,n_plot=0):
# Look through file to find number of trajectories
  raw_data = open(datafile, 'r')
  traj_no = 1
  for line in raw_data:
    if line.strip() == '':     # blank lines divide different trajectories
      traj_no += 1
  raw_data.close()
  print str(traj_no)+' trajectories detected'

  if spacetime:                # 3d plot requested with time along z axis
    pass                       # Not currently implemented
  else:
    plotstr = 'plot \''+datafile+'\' every :::0::0 with line ls 1 notitle'
    if traj_no > 1:            # more than one trajectory
      for i in range(1,traj_no):
        plotstr += ',\\\n \''+datafile+'\' every :::'+str(i)+'::'+str(i)+' with\
   line lc '+str(i)+' notitle'
    plotstr += '\n'
    return plotstr

def plot_2d_map(datafile,colour,shading,contours,contour_no,linegradient,
n_plot=0):
  plotstr = ''
  if contours:
    try:
      contour_datafile = datafile.split('_')[0]+'_contour_'+datafile.split('_')\
[2]
    except IndexError:
      contour_datafile = 'contour_data_'+str(n_plot)+'.dat'
# Generate file containing contour data
    plotstr += 'set contour\n'
    plotstr += 'set cntrparam levels auto '+str(contour_no)+'\n'
# Output contour data to file
    plotstr += 'unset surface\n'
    plotstr += 'set table \''+contour_datafile+'\'\n'
    plotstr += 'splot '+'\''+datafile+'\' using 1:2:4 with pm3d\n'
    plotstr += 'unset table\n'
  if shading and contours:           # contours on top of a shaded map
    plotstr += 'plot '+'\''+datafile+'\' using 1:2:4 with image notitle, \''+\
contour_datafile+'\' with line ls 1 notitle\n'
  elif shading:                      # shading only
    plotstr += 'plot '+'\''+datafile+'\' using 1:2:4 with image notitle\n'
  elif contours and linegradient:    # contours with coloured lines
    plotstr += 'plot '+'\''+contour_datafile+'\' with line lc palette notitle\n'
  else:                              # contours with black lines
    plotstr += 'plot '+'\''+contour_datafile+'\' with line ls 1 notitle\n'
  return plotstr


def plot_2d_density(datafile,colour,shading,mesh,linegradient,n_plot=0):
  plotstr = ''
  if shading and mesh:        # mesh on top of a shaded surface
    plotstr += 'set multiplot\n'
    plotstr += 'unset hidden3d\n'
    plotstr += 'splot \''+datafile+'\' using 1:2:4 with pm3d notitle\n'
    plotstr += 'set xrange restore\n'
    plotstr += 'set yrange restore\n'
    plotstr += 'set zrange restore\n'
    plotstr += 'set cbrange restore\n'
    plotstr += 'unset ztics\n'
    plotstr += 'unset xlabel\n'
    plotstr += 'unset ylabel\n'
    plotstr += 'set hidden3d\n'
    plotstr += 'splot \''+datafile+'\' using 1:2:4 every ev:ev with\
 line lc \'white\' lw 0.5 notitle\n'
    plotstr += 'unset multiplot\n'
  elif mesh:                   # mesh only
    plotstr += 'set hidden3d\n'
    if linegradient:           # - coloured lines
      plotstr += 'splot '+'\''+datafile+'\' using 1:2:4 every ev:ev\
 with line lc palette lw  0.5 notitle\n'
    else:                      # - black lines
      plotstr += 'splot '+'\''+datafile+'\' using 1:2:4 every ev:ev\
 with line ls 1 lw 0.5 notitle\n'
  else:                        # shading only
    plotstr += 'splot '+'\''+datafile+'\' using 1:2:4 with pm3d notitle\n'
    
  return plotstr


def plot_2d_bar(datafile,colour,shading,mesh,linegradient,n_plot=0):
  plotstr = ''
  try:
    bar_datafile = datafile.split('_')[0]+'_bar_'+datafile.split('_')[2]
  except IndexError:
    bar_datafile = 'bar_data_'+str(n_plot)+'.dat'
# Generate file containing bar data
  raw_datafile = open(datafile, 'r')
  line = raw_datafile.readline()
  while line.strip()[0] == '#':   # read past comments
    line = raw_datafile.readline()
  point = [float(x) for x in line.split()]
  x1 = point[0]                   # calculate point spacing
  y1 = point[1]
  point = [float(x) for x in raw_datafile.readline().split()]
  x2 = point[0]
  while line.strip() != '':
    line = raw_datafile.readline()
  point = [float(x) for x in raw_datafile.readline().split()]
  y2 = point[1]
  delta_x = 0.5 * (x2 - x1)
  delta_y = 0.5 * (y2 - y1)
  raw_datafile.seek(0)            # generate plot file with bar data
  plotfile = open(bar_datafile, 'w')
  block_ym = []
  block_yp = []
  for line in raw_datafile:
    if line[0] == '#':
      pass
    elif line.strip() == '':
      for point in block_ym:
        plotfile.write(str(point[0])+' '+str(point[1])+' '+str(point[2])+'\n')
      plotfile.write('\n')
      block_ym = []
      for point in block_yp:
        plotfile.write(str(point[0])+' '+str(point[1])+' '+str(point[2])+'\n')
      plotfile.write('\n')
      block_yp = []
    else:
      point = [float(x) for x in line.split()]
      block_ym.append([point[0]-delta_x, point[1]-delta_y, point[3]])
      block_ym.append([point[0]+delta_x, point[1]-delta_y, point[3]])
      block_yp.append([point[0]-delta_x, point[1]+delta_y, point[3]])
      block_yp.append([point[0]+delta_x, point[1]+delta_y, point[3]])
  raw_datafile.close()
  plotfile.close()

  if mesh and shading:
    plotstr += 'set multiplot\n'
    plotstr += 'unset hidden3d\n'
    plotstr += 'splot \''+bar_datafile+'\' with pm3d notitle\n'
    plotstr += 'set hidden3d\n'
    plotstr += 'set xrange restore\n'
    plotstr += 'set yrange restore\n'
    plotstr += 'set zrange restore\n'
    plotstr += 'set cbrange restore\n'
    plotstr += 'unset ztics\n'
    plotstr += 'unset xlabel\n'
    plotstr += 'unset ylabel\n'
    plotstr += 'splot \''+bar_datafile+'\' with lines ls 1 lw 0.5 notitle\n'
    plotstr += 'unset multiplot\n'
  elif mesh:
    plotstr += 'set hidden3d\n'
    if linegradient:
      plotstr += 'splot \''+bar_datafile+'\' with lines lc palette lw 0.5\
 notitle\n'
    else:
      plotstr += 'splot \''+bar_datafile+'\' with lines ls 1 lw 0.5 notitle\n'
  else:
      plotstr += 'splot \''+bar_datafile+'\' with pm3d\n'
  return plotstr


# ----------------------------- 3D PLOT OPTIONS ----------------------------- #

def plot_3d_trajectory(datafile,n_plot=0):
  plotstr = ''
  # Look through file to find number of trajectories
  raw_data = open(datafile, 'r')
  traj_no = 1
  for line in raw_data:
    if line.strip() == '':     # blank lines divide different trajectories
      traj_no += 1
  raw_data.close()
  print str(traj_no)+' trajectories detected'

  plotstr = 'splot \''+datafile+'\' every :::0::0 with line ls 1 notitle'
  if traj_no > 1:              # more than one trajectory
    for i in range(1,traj_no):
      plotstr += ',\\\n \''+datafile+'\' every :::'+str(i)+'::'+str(i)+' with\
 line ls '+str(i)+' notitle'
  plotstr += '\n'
  return plotstr

def plot_3d_sweep(datafile,plane,no3d,gpif,n_plot=0):
  from plot_utils import settitle
  opt = gpif.get('plot_title')
  autoscale = gpif.get('autoscale_'+plane,'boolean')
  if autoscale:
    pass
  else:
    amin = gpif.get(plane+'min','float')
    amax = gpif.get(plane+'max','float')
  title = settitle(opt,datafile)
  plotstr = ''
  raw_data = open(datafile,'r')
  line = raw_data.readline().strip()
  n_lattice = 0
  lat_point = []
  while line != '':
    if line[0] == '#':
      pass
    else:
      n_lattice += 1
      lat_point.append(float(line.split()[0]))
    line = raw_data.readline().strip()
  start = 0
  stop = n_lattice
  if not autoscale:
    i = 1
    for point in lat_point:
      if point < amin:
        start = i 
      elif point < amax:
        stop = i
      i += 1
  for i in range(start,stop):
    if title[0] == 'u':
      pass
    else:
      plotstr += title[:-2]+', '+plane+' = '+str(lat_point[i])[:6]+'\'\n'
    if no3d:
      plotstr += 'set size ratio -1\nplot \''+datafile+'\' '
      if plane == 'x':
        plotstr += 'using 2:3:4 every ::'+str(i)+'::'+str(i)+' with\
 image notitle\n'
      elif plane == 'y':
        plotstr += 'using 1:3:4 every :32::'+str(i)+' with image\
 notitle\n'
      elif plane == 'z':
        plotstr += 'using 1:2:4 every :::'+str(n_lattice*i)+'::'+\
  str(n_lattice*(i+1)-1)+' with image notitle\'\n'
      else:
        print 'Invalid plane: '+plane+' for 3D sweep plot'
        exit()
    else:
      plotstr += 'splot \''+datafile+'\' using 1:2:3:4 '
      if plane == 'x':
        plotstr += 'every ::'+str(i)+'::'+str(i)+' with image notitle\n'
      elif plane == 'y':
        plotstr += 'every :32::'+str(i)+' with image notitle\n'
      elif plane == 'z':
        plotstr += 'every :::'+str(n_lattice*i)+'::'+str(n_lattice*(i+1)-1)+\
  ' with image notitle\'\n'
      else:
        print 'Invalid plane: '+plane+' for 3D sweep plot'
        exit()  
  return plotstr


