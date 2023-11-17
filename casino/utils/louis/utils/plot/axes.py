#!/usr/bin/python
#------------------------------------------------------------------------------#
# Functions to set axes for the different types of plot available. All return  #
# a string containing the gnuplot commands for setting the relevant axes       #
#------------------------------------------------------------------------------#

import sys
import math
from plot_utils import cont

# --------------------------------- SYNC ------------------------------------ #

def syncaxes():
  print 'Autoscaling requested for multiple datafiles. Apply the same scaling\
 to all\nplots? (y\\n)'
  if cont():
    return True
  else:
    return False

# ---------------------------- AUTOSCALE FUNCTION --------------------------- # 

def autoscale(datafiles,indes,deps,conts):
  ranges = 2*(len(indes)+len(deps))*[0]
  n_file = 0
  for datafile in datafiles:
    raw_data = open(datafile,'r')
    line = raw_data.readline()
    while line[0] == '#':
      line = raw_data.readline()
    point = [float(x) for x in line.strip().split()]
    if n_file == 0:
      i = 0
      for col in indes:
        ranges[i] = point[col]
        ranges[i+1] = point[col]
        i += 2
      for line in raw_data:
        i = 0
        point = [float(x) for x in line.strip().split()]
        if not point:
          continue
        for col in indes:
          if point[col] < ranges[i]:
            ranges[i] = point[col]
          elif point[col] > ranges[i+1]:
            ranges[i+1] = point[col]
          i += 2
        usepoint = True
        for axis in conts:
          if point[axis[0]] < axis[1] or point[axis[0]] > axis[2]:
            usepoint = False
        if usepoint:
          for col in deps:
            ranges[i] = point[col]
            ranges[i+1] = point[col]
            i += 2
          break
    for line in raw_data:
      i = 0
      point = [float(x) for x in line.strip().split()]
      if not point:
        continue
      for col in indes:
        if point[col] < ranges[i]:
          ranges[i] = point[col]
        elif point[col] > ranges[i+1]:
          ranges[i+1] = point[col]
        i += 2
      usepoint = True
      for axis in conts:
        if not (point[axis[0]] > axis[1] and point[axis[0]] < axis[2]):
          usepoint = False  
      if usepoint:
        for col in deps:
          if point[col] < ranges[i]:
            ranges[i] = point[col]
          elif point[col] > ranges[i+1]:
            ranges[i+1] = point[col]
          i += 2
    n_file += 1
  for j in range(len(indes)+len(deps)):
    i = 2*j
    try:
      alrange = math.log((ranges[i+1]-ranges[i]),2)
    except ValueError:
      print i,(ranges[i+1]-ranges[i])
      alrange = 0
    a_inter = 2**(math.floor(alrange-2.585))
    ranges[i] = a_inter * math.floor(ranges[i]/a_inter)
    ranges[i+1] = a_inter * math.ceil(ranges[i+1]/a_inter)
  return ranges

# ----------------------------- 1D AXIS OPTIONS ----------------------------- #

def setaxes_1d_trajectory(gpif,datafiles,rec=False,sync=False):
  axisstr = ''
  if rec:
    return axisstr
  axislabels = gpif.get('axislabels','boolean')
  auto_t = gpif.get('autoscale_t','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  ttic = gpif.get('tticks','float')
  xtic = gpif.get('xticks','float')
  if not auto_t:
    tmin = gpif.get('tmin','float')
    tmax = gpif.get('tmax','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  autosc = auto_x or auto_t
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    if not auto_t:
      axisstr += 'set xrange ['+str(tmin)+':'+str(tmax)+']\n'
    if not auto_x:
      axisstr += 'set yrange ['+str(xmin)+':'+str(xmax)+']\n'
  else:        # Same scaling for all plots
    if auto_t and auto_x:
      tmin,tmax,xmin,xmax = autoscale(datafiles,[0,1],[],[])
      print 't range autoscaled to ['+str(tmin)+':'+str(tmax)+']'
      print 'x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
    elif auto_t:
      tmin,tmax = autoscale(datafiles,[0],[],[])
      print 't range autoscaled to ['+str(tmin)+':'+str(tmax)+']'
    elif auto_x:
      xmin,xmax = autoscale(datafiles,[],[1],[(0,tmin,tmax)])
      print 'x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
    axisstr += 'set xrange ['+str(tmin)+':'+str(tmax)+']\n'
    axisstr += 'set yrange ['+str(xmin)+':'+str(xmax)+']\n'
  if ttic == 0:
    axisstr += 'unset xtics\n'
  elif ttic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(ttic)+'\n'
  if xtic == 0:
    axisstr += 'unset ytics\n'
  elif xtic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'t\'\n'
    axisstr += 'set ylabel \'x\'\n'
  return (sync,axisstr)

def setaxes_1d_density(gpif,datafiles,rec=False,sync=False):
  axisstr = ''
  if rec:
    return axisstr
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  autosc = auto_x or auto_y
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
  else:        # Same scaling for all plots
    if auto_y and auto_x:
      xmin,xmax,ymin,ymax = autoscale(datafiles,[0,1],[],[])
      print 'x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
      print 'y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    elif auto_x:
      xmin,xmax = autoscale(datafiles,[0],[],[])
      print 'x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
    elif auto_y:
      ymin,ymax = autoscale(datafiles,[],[1],[(0,xmin,xmax)])
      print 'y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
  return (sync,axisstr)

# ----------------------------- 2D AXIS OPTIONS ----------------------------- #

def setaxes_2d_trajectory(gpif,datafiles,rec=False,sync=False):
  if rec:
    return ''
  axisstr = 'set size ratio -1\n' # 1 unit on x axis = 1 unit on y axis
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  autosc = auto_x or auto_y
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
  else:        # Same scaling for all plots
    if auto_y and auto_x:
      xmin,xmax,ymin,ymax = autoscale(datafiles,[0,1],[],[])
      print 'x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
      print 'y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    elif auto_x:
      xmin,xmax = autoscale(datafiles,[],[0],[(1,ymin,ymax)])
      print 'x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
    elif auto_y:
      ymin,ymax = autoscale(datafiles,[],[1],[(0,xmin,xmax)])
      print 'y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
    axisstr += 'set ylabel \'y\'\n'
  return (sync,axisstr)

def setaxes_2d_spacetime(gpif,datafiles,rec=False,sync=False):
  pass

def setaxes_2d_map(gpif,datafiles,rec=False,sync=False):
  if rec:
    return ''
  axisstr = 'set size ratio -1\n' # 1 unit on x axis = 1 unit on y axis
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  auto_c = gpif.get('autoscale_z','boolean')
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  ctic = gpif.get('zticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  if not auto_c:
    cmin = gpif.get('zmin','float')
    cmax = gpif.get('zmax','float')
  autosc = auto_x or auto_y or auto_c
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
    if not auto_c:
      axisstr += 'set cbrange ['+str(cmin)+':'+str(cmax)+']\n'
  else:        # Same scaling for all plots
    if auto_c:
      if auto_y and auto_x:
        xmin,xmax,ymin,ymax,cmin,cmax = autoscale(datafiles,[0,1,3],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
        print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
      elif auto_x:
        xmin,xmax,cmin,cmax = autoscale(datafiles,[0],[3],[(1,ymin,ymax)])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
      elif auto_y:
        ymin,ymax,cmin,cmax = autoscale(datafiles,[1],[3],[(0,xmin,xmax)])
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
        print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
      else:
        cmin,cmax = autoscale(datafiles,[],[3],[(0,xmin,xmax),(1,ymin,ymax)])
        print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
    else:
      if auto_y and auto_x:
        xmin,xmax,ymin,ymax = autoscale(datafiles,[0,1],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
      elif auto_x:
        xmin,xmax = autsocale(datafiles,[0],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
      elif auto_y:
        ymin,ymax = autoscale(datafiles,[1],[],[])
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
    axisstr += 'set cbrange ['+str(cmin)+':'+str(cmax)+']\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
   axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if ctic == 0:
    axisstr += 'unset cbtics\n'
  elif ctic < 0:
    pass
  else:
    axisstr += 'set cbtics '+str(ctic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
    axisstr += 'set ylabel \'y\'\n'
  return (sync,axisstr)

def setaxes_2d_density(gpif,datafiles,rec=False,sync=False):
  axisstr = 'set size ratio -1\n' # 1 unit on x axis = 1 unit on y axis
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  auto_z = gpif.get('autoscale_z','boolean')
  if rec:
    if not sync:
      axisstr = ''
      if auto_x:
        axisstr += 'set xrange [*:*] writeback\n'
      if auto_y:
        axisstr += 'set yrange [*:*] writeback\n'
      if auto_z:
        axisstr += 'set zrange [*:*] writeback\n'
      axisstr += 'set cbrange [*:*] writeback\n'
    mesh = gpif.get('mesh','boolean')
    shading = gpif.get('shading','boolean')
    plot_type = gpif.get('plot_type')
    if (plot_type == 'density') and mesh:
      if not auto_x:
        xmin = gpif.get('xmin','float')
        xmax = gpif.get('xmax','float')
      if not auto_y:
        ymin = gpif.get('ymin','float')
        ymax = gpif.get('ymax','float')
      nx = 0
      ny = 0
      raw_data = open(datafiles[0],'r')
      line = raw_data.readline()
      while line[0] == '#':
        line = raw_data.readline()
      point = [float(x) for x in line.split()]
      if auto_y:
        ny += 1
      elif point[1] > ymin and point[1] < ymax:
        ny += 1
      while line:
        if auto_x:
          nx += 1
        elif point[0] > xmin and point[0] < xmax:
          nx += 1
        line = raw_data.readline().strip()
        point = [float(x) for x in line.split()]
      line = raw_data.readline()
      point = [float(x) for x in line.split()]
      if auto_y:
        ny += 1
      elif point[1] > ymin and point[1] < ymax:
        ny += 1
      while line:
        if not line.strip():
          line = raw_data.readline()
          yval = float(line.split()[1])
          if auto_y:
            ny += 1
          elif yval > ymin and yval < ymax:
            ny += 1
        else:
          line = raw_data.readline()
      ev = int(math.ceil(max(float(nx)/40.0,float(ny)/40.0)))
      axisstr += 'ev = '+str(ev)+'\n'
      if mesh and shading:
        ztic = gpif.get('zticks','float')
        if ztic > 0:
          axisstr += 'set ztics out '+str(ztic)+'\n'
        elif ztic == 0:
          axisstr += 'unset ztics\n'
        else:
          axisstr += 'set ztics out\n'
        if axislabels:
          axisstr += 'set xlabel \'x\'\n'
          axisstr += 'set ylabel \'y\'\n'
    return axisstr
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  ztic = gpif.get('zticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  if not auto_z:
    zmin = gpif.get('zmin','float')
    zmax = gpif.get('zmax','float')
  autosc = auto_x or auto_y or auto_z
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+'] writeback\n'
    else:
      axisstr += 'set xrange [] writeback\n'
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+'] writeback\n'
    else:
      axisstr += 'set yrange [] writeback\n'
    if not auto_z:
      axisstr += 'set zrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
      axisstr += 'set cbrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
    else:
      axisstr += 'set zrange [] writeback\n'
      axisstr += 'set cbrange [] writeback\n'
  else:        # Same scaling for all plots
    if auto_z:
      if auto_y and auto_x:
        xmin,xmax,ymin,ymax,zmin,zmax = autoscale(datafiles,[0,1,3],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
        print 'cb range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
      elif auto_x:
        xmin,xmax,zmin,zmax = autoscale(datafiles,[0],[3],[(1,ymin,ymax)])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
        print 'cb range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
      elif auto_y:
        ymin,ymax,zmin,zmax = autoscale(datafiles,[1],[3],[(0,xmin,xmax)])
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
        print 'cb range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
      else:
        zmin,zmax = autoscale(datafiles,[],[3],[(0,xmin,xmax),(1,ymin,ymax)])
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
        print 'cb range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
    else:
      if auto_y and auto_x:
        xmin,xmax,ymin,ymax = autoscale(datafiles,[0,1],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
      elif auto_x:
        xmin,xmax = autsocale(datafiles,[0],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
      elif auto_y:
        ymin,ymax = autoscale(datafiles,[1],[],[])
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+'] writeback\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+'] writeback\n'
    axisstr += 'set zrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
    axisstr += 'set cbrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if ztic == 0:
    axisstr += 'unset ztics\n'
  elif ztic < 0:
    axisstr += 'set ztics out autofreq\n'
  else:
    axisstr += 'set ztics out '+str(ztic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
    axisstr += 'set ylabel \'y\'\n'
  return (sync,axisstr)

# ----------------------------- 3D AXIS OPTIONS ----------------------------- #

def setaxes_3d_trajectory(gpif,datafiles,rec=False,sync=False):
  axisstr = ''
  if rec:
    return axisstr
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  auto_z = gpif.get('autoscale_z','boolean')
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  ztic = gpif.get('zticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  if not auto_z:
    zmin = gpif.get('zmin','float')
    zmax = gpif.get('zmax','float')
  autosc = auto_x or auto_y or auto_z
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+'] writeback\n'
    else:
      axisstr += 'set xrange [] writeback\n'
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+'] writeback\n'
    else:
      axisstr += 'set yrange [] writeback\n'
    if not auto_z:
      axisstr += 'set zrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
      axisstr += 'set cbrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
    else:
      axisstr += 'set zrange [] writeback\n'
      axisstr += 'set cbrange [] writeback\n'
  else:        # Same scaling for all plots
    if auto_z:
      if auto_x and auto_y:
        xmin,xmax,ymin,ymax,zmin,zmax = autoscale(datafiles,[0,1,2],[],[])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
      elif auto_x:
        xmin,xmax,zmin,zmax = autoscale(datafiles,[0],[2],[(1,ymin,ymax)])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
      elif auto_y:
        ymin,ymax,zmin,zmax = autoscale(datafiles,[1],[2],[(0,xmin,xmax)])
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
      else:
        zmin,zmax = autoscale(datafiles,[],[2],[(0,xmin,xmax),(1,ymin,ymax)])
        print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
    else:
      if auto_y and auto_x:
        xmin,xmax,ymin,ymax = autoscale(datafiles,[],[0,1],[(2,zmin,zmax)])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
      elif auto_x:
        xmin,xmax = autsocale(datafiles,[],[0],[(1,ymin,ymax),(2,zmin,zmax)])
        print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
      elif auto_y:
        ymin,ymax = autoscale(datafiles,[],[1],[(0,xmin,xmax),(2,zmin,zmax)])
        print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+'] writeback\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+'] writeback\n'
    axisstr += 'set zrange ['+str(zmin)+':'+str(zmax)+'] writeback\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if ztic == 0:
    axisstr += 'unset ztics\n'
  elif ztic < 0:
    axisstr += 'set ztics out autofreq\n'
  else:
    axisstr += 'set ztics out '+str(ztic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
    axisstr += 'set ylabel \'y\'\n'
    axisstr += 'set zlabel \'z\'\n'
  return (sync,axisstr)

def setaxes_3d_sweep(gpif,datafiles,rec=False,sync=False):
  axisstr = ''
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  auto_z = gpif.get('autoscale_z','boolean')
  auto_c = True
  if rec:
    if sync:
      pass
    else:
      conts = []
      if not auto_x:
        conts.append((0,gpif.get('xmin','float'),gpif.get('xmax','float')))
      if not auto_y:
        conts.append((1,gpif.get('ymin','float'),gpif.get('ymax','float')))
      if not auto_z:
        conts.append((2,gpif.get('zmin','float'),gpif.get('zmax','float')))
      cmin,cmax = autoscale([datafiles[0]],[],[3],conts)
      axisstr += 'set cbrange ['+str(cmin)+':'+str(cmax)+']\n'
      print 'cb range scaled to ['+str(cmin)+':'+str(cmax)+']'
    return axisstr
  axisstr += 'set size ratio -1\n'
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  ztic = gpif.get('zticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  if not auto_z:
    zmin = gpif.get('zmin','float')
    zmax = gpif.get('zmax','float')
  autosc = auto_x or auto_y or auto_z or auto_c # Currently always True...
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    conts = []
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
      conts.append((0,xmin,xmax))
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
      conts.append((1,ymin,ymax))
    if not auto_z:
      conts.append((2,zmin,zmax))
    if not auto_c:
      axisstr += 'set cbrange [+str(cmin)+:+str(cmax)+]\n'
  else:        # Same scaling for all plots
    if auto_c:
      if auto_z:
        if auto_x and auto_y:
          xmin,xmax,ymin,ymax,zmin,zmax,cmin,cmax = autoscale(datafiles,\
[0,1,2,3],[],[])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_x:
          xmin,xmax,zmin,zmax,cmin,cmax = autoscale(datafiles,[0,2],[3],\
[(1,ymin,ymax)])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_y:
          ymin,ymax,zmin,zmax,cmin,cmax = autoscale(datafiles,[1,2],[3],\
[(0,xmin,xmax)])
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        else:
          zmin,zmax,cmin,cmax = autoscale(datafiles,[2],[3],[(0,xmin,xmax),\
(1,ymin,ymax)])
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
      else:
        if auto_x and auto_y:
          xmin,xmax,ymin,ymax,cmin,cmax = autoscale(datafiles,[0,1,3],[],[])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_x:
          xmin,xmax,cmin,cmax = autoscale(datafiles,[0],[3],[(1,ymin,ymax),\
(2,zmin,zmax)])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_y:
          ymin,ymax,cmin,cmax = autoscale(datafiles,[1],[3],[(0,xmin,xmax),\
(2,zmin,zmax)])
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        else:
          cmin,cmax = autoscale(datafiles,[],[3],[(0,xmin,xmax),(1,ymin,ymax),\
(2,zmin,zmax)])
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
    else:
      print 'cb axis not autoscaled... Should not happen'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
    axisstr += 'set cbrange ['+str(cmin)+':'+str(cmax)+']\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
    axisstr += 'set ylabel \'y\'\n'
  return (sync,axisstr)

def setaxes_3d_sweep3d(gpif,datafiles,rec=False,sync=False):
  axisstr = ''
  axislabels = gpif.get('axislabels','boolean')
  auto_x = gpif.get('autoscale_x','boolean')
  auto_y = gpif.get('autoscale_y','boolean')
  auto_z = gpif.get('autoscale_z','boolean')
  auto_c = True
  if rec:
    if sync:
      pass
    else:
      conts = []
      if not auto_x:
        conts.append((0,gpif.get('xmin','float'),gpif.get('xmax','float')))
      if not auto_y:
        conts.append((1,gpif.get('ymin','float'),gpif.get('ymax','float')))
      if not auto_z:
        conts.append((2,gpif.get('zmin','float'),gpif.get('zmax','float')))
      cmin,cmax = autoscale([datafiles[0]],[],[3],conts)
      axisstr += 'set cbrange ['+str(cmin)+':'+str(cmax)+']\n'
      print 'cb range scaled to ['+str(cmin)+':'+str(cmax)+']'
    return axisstr
  xtic = gpif.get('xticks','float')
  ytic = gpif.get('yticks','float')
  ztic = gpif.get('zticks','float')
  if not auto_x:
    xmin = gpif.get('xmin','float')
    xmax = gpif.get('xmax','float')
  if not auto_y:
    ymin = gpif.get('ymin','float')
    ymax = gpif.get('ymax','float')
  if not auto_z:
    zmin = gpif.get('zmin','float')
    zmax = gpif.get('zmax','float')
  autosc = auto_x or auto_y or auto_z or auto_c # Currently always True...
  if autosc and len(datafiles) > 1:
    sync = syncaxes()
  else:
    sync = False
  if not sync: # Scale differently for each plot.
    conts = []
    if not auto_x:
      axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
      conts.append((0,xmin,xmax))
    if not auto_y:
      axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
      conts.append((1,ymin,ymax))
    if not auto_z:
      axisstr += 'set zrange ['+str(zmin)+':'+str(zmax)+']\n'
      conts.append((2,zmin,zmax))
    if not auto_c:
      axisstr += 'set cbrange [+str(cmin)+:+str(cmax)+]\n'
  else:        # Same scaling for all plots
    if auto_c:
      if auto_z:
        if auto_x and auto_y:
          xmin,xmax,ymin,ymax,zmin,zmax,cmin,cmax = autoscale(datafiles,\
[0,1,2,3],[],[])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_x:
          xmin,xmax,zmin,zmax,cmin,cmax = autoscale(datafiles,[0,2],[3],\
[(1,ymin,ymax)])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_y:
          ymin,ymax,zmin,zmax,cmin,cmax = autoscale(datafiles,[1,2],[3],\
[(0,xmin,xmax)])
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        else:
          zmin,zmax,cmin,cmax = autoscale(datafiles,[2],[3],[(0,xmin,xmax),\
(1,ymin,ymax)])
          print ' z range autoscaled to ['+str(zmin)+':'+str(zmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
      else:
        if auto_x and auto_y:
          xmin,xmax,ymin,ymax,cmin,cmax = autoscale(datafiles,[0,1,3],[],[])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_x:
          xmin,xmax,cmin,cmax = autoscale(datafiles,[0],[3],[(1,ymin,ymax),\
(2,zmin,zmax)])
          print ' x range autoscaled to ['+str(xmin)+':'+str(xmax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        elif auto_y:
          ymin,ymax,cmin,cmax = autoscale(datafiles,[1],[3],[(0,xmin,xmax),\
(2,zmin,zmax)])
          print ' y range autoscaled to ['+str(ymin)+':'+str(ymax)+']'
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
        else:
          cmin,cmax = autoscale(datafiles,[],[3],[(0,xmin,xmax),(1,ymin,ymax),\
(2,zmin,zmax)])
          print 'cb range autoscaled to ['+str(cmin)+':'+str(cmax)+']'
    else:
      print 'cb axis not autoscaled... Should not happen'
    axisstr += 'set xrange ['+str(xmin)+':'+str(xmax)+']\n'
    axisstr += 'set yrange ['+str(ymin)+':'+str(ymax)+']\n'
    axisstr += 'set zrange ['+str(zmin)+':'+str(zmax)+']\n'
    axisstr += 'set cbrange ['+str(cmin)+':'+str(cmax)+']\n'
  if xtic == 0:
    axisstr += 'unset xtics\n'
  elif xtic < 0:
    axisstr += 'set xtics out nomirror\n'
  else:
    axisstr += 'set xtics out nomirror '+str(xtic)+'\n'
  if ytic == 0:
    axisstr += 'unset ytics\n'
  elif ytic < 0:
    axisstr += 'set ytics out nomirror\n'
  else:
    axisstr += 'set ytics out nomirror '+str(xtic)+'\n'
  if ztic == 0:
    axisstr += 'unset ztics\n'
  elif ztic < 0:
    axisstr += 'set ztics out autofreq\n'
  else:
    axisstr += 'set ztics out '+str(ztic)+'\n'
  if axislabels:
    axisstr += 'set xlabel \'x\'\n'
    axisstr += 'set ylabel \'y\'\n'
    axisstr += 'set zlabel \'z\'\n'
  return (sync,axisstr)



