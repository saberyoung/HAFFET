#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : src/sneGui.py
# Author            : syang <sheng.yang@astro.su.se>
# Date              : 12.11.2019
# Last Modified Date: 11.02.2022
# Last Modified By  : syang <sheng.yang@astro.su.se>

from __future__ import print_function
from builtins import input
import os, sys, glob, warnings, time
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from tabulate import tabulate
from weakref import ref
from astropy.time import Time

from sdapy import __version__, snerun, \
   image_tool, constants, read_default
from sdapy.pbar import get_progress_bar
from sdapy.functions import *
from sdapy.filters import *
from sdapy.model_fitters import get_pars, get_model, get_engine
from sdapy.models.arnett_tail import functions as radioactivemodels
from sdapy.models.sbo import functions as sbomodels
from sdapy.models.bazin.functions import bazin as bazinfunc
from sdapy.models.risepl.functions import powerlaw_full as plfunc

#from distutils.spawn import find_executable
import matplotlib.pyplot as plt
# Check if latex is installed
#if find_executable('latex'):
#   plt.rcParams['text.usetex'] = True
#else:
# disable latex in matplotlib
# if latex working properly, set it to True
plt.rcParams["text.usetex"] = False

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Cursor
from matplotlib.transforms import CompositeGenericTransform
from matplotlib.gridspec import SubplotSpec

# tkinter for the gui
if sys.version_info.major == 3: 
   import tkinter as Tk
   from tkinter import filedialog, font
else:
   import Tkinter as Tk
   import tkFileDialog as filedialog
   import tkFont as font
from tkinter import messagebox

dpi = 96
dataoptions = {   
   'Data preparation': 'down/re-load ZTF/ATLAS/externel forced/alert photometry/spectra',
   'Multiband LC fits': 'Fit multiband LCs, e.g. via Bazin fits for main peak, power law for early LCs, etc, with Gaussian Process (george), Monte Carlo (emcee), or minimize (scipy)',
   'SED fits': 'calculate magnitude differences bwtween different bands/phases.\n\nestimate host ebv by comparing to intrinstic colors\n\nset peak/explosion epochs.\n\ncalculate bolometric LCs with Blackbody fits, or Lyman 2016 analytic functions',
   'Bolometric LC fits': 'Fit bolometric LCs, e.g. shock cooling breakout fits for early tail, Arnett fits for the main peak, tail fits, etc.',
   'Spectral line fits': 'Fit for spectral element line emission/absorption, e.g. via Gaussian or viglot',
}

plotoptions = {
   'Summarize Plot': 'Summarize plot',   
   'Flux LC': 'Multi-band Lightcurves in uJy',
   'Mag LC' : 'Multi-band Lightcurves in AB magnitude',   
   'colour': 'Colour evlotion plots',
   'luminosity' : 'Bolometric Lightcurves in erg/s or magnitude',
   'spectra' : 'Spectral full range or specific line range',   
   'separator1': None,
   'Finder': 'create a finder with online sources, e.g. PS1',
   'LC baseline': 'Baseline check plots for ZTF forced photometry',         
   'SNe scatter': 'Population for 1 parameter (Histogram) or 2 parameters (scatter), show their correlation',   
}

# ztfquery source path
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

###############################     MAIN WINDOW ##########################
class MainMenu(Tk.Frame):

   def __init__(self, parent, width, height):      

      # check dir
      if not snerun.check_dir(check=True): sys.exit ('check data dir first')
      
      global snelist_pars, snobject_pars, searchnames
      snelist_pars = read_default.get_parameters(keylist='snelist')['snelist']   
      snobject_pars = read_default.get_parameters(keylist='snobject')['snobject']
      if len(snelist_pars['idkey1']) > 0: searchnames = [snelist_pars['idkey1']]
      else:  searchnames = []
      
      Tk.Frame.__init__(self, parent)
      self.parent = parent
      self.width = width
      self.height = height
      
      # init sdapy class      
      self.obj = Tk.StringVar()
      self.obj.set(None)      
      self.ztfp = snerun.snelist()
      
      # initialize
      self.ztfp.parse_meta()
      self.ztfp.parse_params()
      
      # init GUI
      self.initUI()
      self.Message( snerun.print_logo(returns=True) )
      
   def initUI(self):    ##########################################

      """ INIT Frames: divide Parent Menu """
      
      """ ROW 1: SN meta (Left span) """ 
      row1left = Tk.Frame(self.parent)      
      row1left.grid(row=0,column=0,sticky=Tk.N,padx=5,pady=5)

      ### split row1 left into 3 rowspans
      row1left.rowconfigure(0, weight=1) # obj label
      row1left.rowconfigure(1, weight=2) # objlist search & input
      row1left.rowconfigure(2, weight=2) # data
      row1left.rowconfigure(3, weight=2) # plot
      row1left.rowconfigure(4, weight=2) # fit
      row1left.rowconfigure(5, weight=2) # config & quit
      
      ''' obj Label frame '''
      labelFrame = Tk.LabelFrame(row1left, text="")
      labelFrame.grid(row=0,column=0,sticky=Tk.NW)
      labelFrame.columnconfigure(0, weight=3)
      labelFrame.columnconfigure(1, weight=1)
      labelFrame.columnconfigure(2, weight=1)
      self.metainfo = Tk.StringVar()
      self.metainfo.set('Obj: %s \n(%i of %i)' % (self.obj.get(), len(self.ztfp.data), len(self.ztfp.meta)))
      tb = Tk.Label(labelFrame,textvariable=self.metainfo,fg="red")      
      tb.grid(row=0,column=0,rowspan=2)

      # search frame
      sFrame = Tk.LabelFrame(labelFrame,text="Query syntax")
      sFrame.grid(row=3,column=0,rowspan=1)
      
      self.sqlcmd = Tk.Entry(sFrame, justify=Tk.LEFT)
      syntax = snelist_pars['syntax']      
      self.sqlcmd.insert(0, syntax)
      self.sqlcmd.grid(row=0, column=0)
      
      def showmeta():
         # update meta
         syntax = self.sqlcmd.get()
         source = self.metasourcevar.get()
         metafile = None
         if not source in ['BTS', 'OAC']:            
            metafile = source
            source = None
         self.ztfp.parse_meta(withnew='new', source=source, metafile=metafile, syntax=syntax, verbose=False)
         self.metainfo.set('Obj: %s \n(%i of %i)' %
                           (self.obj.get(), len(self.ztfp.data), len(self.ztfp.meta)))
         tb = Tk.Label(labelFrame,textvariable=self.metainfo,fg="red")
         # what to show
         showmeta = self.ztfp.meta
         if len( self.ztfp.meta ) > 100:  showmeta = showmeta.head(100)
         self.Message( tabulate(showmeta, headers = 'keys', tablefmt = 'grid') )
      # meta button
      checkMetaBut = Tk.Button(labelFrame,text="Meta",command=showmeta)
      checkMetaBut.grid(row=0, column=1, sticky=Tk.NW)

      # load cache button
      lcBut = Tk.Button(labelFrame,text="load",\
               fg="green", command=self.loadall, activebackground='slate grey')
      lcBut.grid(row=1, column=1, sticky=Tk.NW)
      
      # save cache button
      scBut = Tk.Button(labelFrame,text="save",\
               fg="gold", command=self.saveall, activebackground='slate grey')
      scBut.grid(row=2, column=1, sticky=Tk.NW)  

      # source dropdown      
      def sets(source):
         if source in ['BTS', 'OAC']:
            self.metasourcevar.set(source)
            self.Message('Use meta data from: %s' % source)
         else:
            filetypes = [
               ('text files ', '*.txt'),
               ("CSV Files", "*.csv"),
            ]
            filename = filedialog.askopenfilename(parent=labelFrame, initialdir=os.getcwd(),
                           title='Open the Metafile', filetypes=filetypes)            
            if filename:
               self.metasourcevar.set(filename)
               self.Message('Use meta data from: %s' % filename)         
      sFrame = Tk.LabelFrame(labelFrame,text="source")
      sFrame.grid(row=3,column=1,sticky=Tk.NW)
      self.metasourcevar = Tk.StringVar()
      self.metasourcevar.set('BTS')
      sMenu = Tk.Menubutton(sFrame,textvariable=self.metasourcevar, activebackground='slate grey', width=4)
      sMenu.grid(row=0,column=0)
      picks = Tk.Menu(sMenu)      
      sMenu.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in ['BTS', 'OAC', 'Other']:
         picks.add_command(label=t, command=lambda k=t: sets(k))
         
      ''' objlist frame'''
      objlistFrame = Tk.LabelFrame(row1left, text="Object List")
      objlistFrame.grid(row=1,column=0,sticky=Tk.NW)
      
      # object list
      objlist = [None]
      self.objlistMenu = Tk.OptionMenu(objlistFrame, self.obj, *objlist)
      self.objlistMenu.config(justify=Tk.LEFT,width=8)
      self.objlistMenu.grid(row=0,column=0, sticky=Tk.SW)
      
      # search button
      searchFrameBut = Tk.Button(objlistFrame,text="Search",command=self.SearchObj)
      searchFrameBut.grid(row=1, column=1, sticky=Tk.SW)
      
      # search txt
      oFrame = Tk.LabelFrame(objlistFrame, text="string")
      oFrame.grid(row=1,column=0,sticky=Tk.SW)
      
      self.searchFrameTxt = Tk.Entry(oFrame, justify=Tk.LEFT, width=12)
      self.searchFrameTxt.grid(row=0, column=0, sticky=Tk.SW)      
      
      # add/del
      inputFrameBut = Tk.Button(objlistFrame,text="Add",command=self.addobj)
      inputFrameBut.grid(row=2, column=1, sticky=Tk.SW)
      inputFrameBut = Tk.Button(objlistFrame,text="Del",command=self.delobj)
      inputFrameBut.grid(row=3, column=1, sticky=Tk.SW)
      
      def delcache():
         objid = self.obj.get()
         if objid in ['','None']:
            self.Message('Warning: No objects left')
            return
         f = '%s/cache/%s_data.clf' % (LOCALSOURCE, objid)
         if os.path.exists(f):         
            os.remove(f)
            self.Message('Had removed: %s'%f)
      inputFrameBut = Tk.Button(objlistFrame,text="Remove",fg="black",command=delcache)
      inputFrameBut.grid(row=4, column=1, sticky=Tk.SW)
      
      def goback():
         objid = self.obj.get()
         if objid in ['','None']:
            self.Message('Warning: set object first')
            return
         objlist = np.array(sorted(self.ztfp.data.keys()))
         __ = np.where(objlist == objid)[0][0]         
         if __ == 0:
            self.Message('Warning: %s is already the first object'%objid)
            return
         objid = objlist[__-1]
         self.setobjid(objid)
      inputFrameBut = Tk.Button(objlistFrame,text="\u23EA",command=goback)
      inputFrameBut.grid(row=2, column=0, sticky=Tk.SW)

      def goon():
         objid = self.obj.get()
         if objid in ['','None']:
            self.Message('Warning: set object first')
            return
         objlist = np.array(sorted(self.ztfp.data.keys()))
         __ = np.where(objlist == objid)[0][0]  
         if __ == len(objlist)-1:
            self.Message('Warning: %s is already the last object'%objid)
            return
         objid = objlist[__+1]
         self.setobjid(objid)
      inputFrameBut = Tk.Button(objlistFrame,text="\u23E9",command=goon)
      inputFrameBut.grid(row=2, column=0, sticky=Tk.SE)

      def gorandom():         
         objlist = np.array(sorted(self.ztfp.data.keys()))
         if len(objlist)==0:
            self.Message('Warning: add some objects first')
            return         
         objid = np.random.choice(objlist)          
         self.setobjid(objid)
      inputFrameBut = Tk.Button(objlistFrame,text="\u23ED",command=gorandom)
      inputFrameBut.grid(row=3, column=0, sticky=Tk.SE)  
      
      ''' data frame '''
      def doaction(actstyle): 
         self.optvar.set(actstyle)
         self.Message('Action:\n\n%s: \n%s' % (actstyle, dataoptions[actstyle]))         
      dataFrame = Tk.LabelFrame(row1left,text="For obj")
      dataFrame.grid(row=2,column=0,sticky=Tk.NW)
      self.optvar = Tk.StringVar()
      self.optvar.set(None)
      actionMenu = Tk.Menubutton(dataFrame,textvariable=self.optvar, activebackground='slate grey', width=7)
      actionMenu.grid(row=0,column=0)
      picks   = Tk.Menu(actionMenu)               
      actionMenu.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in dataoptions.keys():
         if dataoptions[t] is None:
            picks.add_separator()
         else:
            picks.add_command(label=t, command=lambda k=t: doaction(k))
            
      ### submit button
      actionButton = Tk.Button(dataFrame,text="submit",command=self.runtasks)
      actionButton.grid(row=0, column=1)
      
      ''' plot frame '''
      plotFrame = Tk.LabelFrame(row1left,text="Plot")
      plotFrame.grid(row=3,column=0,sticky=Tk.NW)
      
      # plot list
      def setplot(figstyle):
         self.Message('show %s plot\n\n%s' % (figstyle, plotoptions[figstyle]))
         self.plotvar.set(figstyle)      
      self.plotvar = Tk.StringVar()
      self.plotvar.set(None)
      plotMenu = Tk.Menubutton(plotFrame,textvariable=self.plotvar, activebackground='slate grey', width=7)
      plotMenu.grid(row=0,column=0)
      picks   = Tk.Menu(plotMenu)
      plotMenu.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in plotoptions.keys():
         if plotoptions[t] is None:
            picks.add_separator()
         else:
            picks.add_command(label=t, command=lambda k=t: setplot(k))
            
      ### submit button
      plotButton = Tk.Button(plotFrame,text="submit",command=self.runplots)
      plotButton.grid(row=0, column=1)
      
      ''' config frame '''      
      configFrame = Tk.LabelFrame(row1left,text="Settings")      
      configFrame.grid(row=4,column=0,sticky=Tk.NW)      
      
      ### clobber and verbose
      self.clobberVar = Tk.IntVar(value=0)
      clobberButton = Tk.Checkbutton(configFrame, text="clobber",variable=self.clobberVar)      
      clobberButton.grid(row=0,column=0,sticky=Tk.W)
      
      self.fitVar = Tk.IntVar(value=1)
      fitButton = Tk.Checkbutton(configFrame, text="Popup",variable=self.fitVar)      
      fitButton.grid(row=0,column=1,sticky=Tk.W)
      
      self.markerVar = Tk.IntVar(value=0)
      markerButton = Tk.Checkbutton(configFrame, text="marker",variable=self.markerVar)
      markerButton.grid(row=1,column=0,sticky=Tk.W)

      self.zoomVar = Tk.IntVar(value=0)
      zoomButton = Tk.Checkbutton(configFrame, text="zoom",variable=self.zoomVar)
      zoomButton.grid(row=1,column=1,sticky=Tk.W)

      # show par
      def showconfig():
         msg = 'for general: \n'
         for key, value in snelist_pars.items(): msg += '%s   :    %s\n' % (key, value)
         msg += '\n\n\nfor single: \n'
         for key, value in snobject_pars.items(): msg += '%s   :    %s\n' % (key, value)
         objid = self.obj.get()
         if objid is None: pass            
         elif objid not in self.ztfp.data: pass
         else:
            msg = 'for general: \n'
            for key, value in snelist_pars.items(): msg += '%s   :    %s\n' % (key, value)
            msg += '\n\n\nfor %s: \n'%objid
            for key, value in self.ztfp.data[objid].kwargs.items(): msg += '%s   :    %s\n' % (key, value)         
         self.Message( msg )
      defBut = Tk.Button(configFrame,text="pars",\
                     fg="blue", command=showconfig, activebackground='slate grey')      
      defBut.grid(row=2, column=0, sticky=Tk.SW)
      
      # quit
      quitBut = Tk.Button(configFrame,text="quit",\
               fg="red", command=self.quit, activebackground='slate grey')      
      quitBut.grid(row=2, column=1, sticky=Tk.SW)
      
      ''' ROW 1: Figure windoe (Right span) '''
      self.row1right = Tk.Frame(self.parent)
      self.row1right.grid(row=0,column=1)      
      
      ### add canvas frame
      self.currentfigure = None
      self.currentfiguresize = None
      cframe = Tk.Frame(self.row1right, name='canvas', bg='white')
      self.fig = Figure(figsize=(self.width/dpi/2, self.height/dpi),
                        tight_layout=True, dpi=dpi)
      self.fig.set_facecolor('white')
      
      self.canvas = FigureCanvasTkAgg(self.fig, master=cframe)
      self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=True)      
      cframe.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)
      
      self.canvas.mpl_connect('scroll_event', self.zoom)
      self.canvas.mpl_connect('button_press_event', self.move_from)      
      self.canvas.mpl_connect('button_release_event', self.move_to)      
      self.canvas.mpl_connect('motion_notify_event', self.motion)
      
      ''' ROW 2: Message menu '''      
      self.row2 = Tk.Frame(self.parent)
      self.row2.grid(row=1,column=0,columnspan=2)
      
      ###  message windows
      self.msgwindow = Tk.Scrollbar(self.row2)
      self.msgwindow.pack(side='right', fill='y')
      self.textbox = Tk.Text(self.row2,yscrollcommand=self.msgwindow.set)
      self.textbox.config(foreground = "red",
            width=int(self.width*.8), height=int(self.height*.02))
      self.textbox.insert(Tk.END,6*"\n"+"!!! WATCH FOR MESSAGES HERE !!!")
      self.msgwindow.config(command=self.textbox.yview)
      self.textbox.config(state='disabled')
      self.textbox.pack()

   ######################################################
   #  functions
   ######################################################
   def set_title(self, ax, title, x=.6, y=1.05):
      
      bbox = dict(boxstyle='square', ec='white', fc='white', alpha=.7)
      ax.text(x,y,title,fontsize=10,weight='bold',transform=ax.transAxes,va='top',ha='right',bbox=bbox)
      
   def set_limits(self, ax, x0=None, xx=None, y0=None, yy=None):
      
      _x0,_xx,_y0,_yy = np.inf,-np.inf,np.inf,-np.inf
      
      for artist in ax.lines + ax.patches:

         if not artist.get_visible(): continue
         if not type(artist.get_transform()) is CompositeGenericTransform: continue

         if hasattr(artist,'get_data'): x,y = map(lambda v:np.append([],v),artist.get_data())
         elif hasattr(artist,'get_xy'):
            try: x,y = map(lambda v:np.append([],v),artist.get_xy().T)
            except:
               self.canvas.draw()
               return         
         if not (len(x) and len(y)): continue
         i = np.arange(len(x))
         
         if not x0 is None: i = i[np.where(x0<=x[i])]
         if not xx is None: i = i[np.where(x[i]<=xx)]
         if not y0 is None: i = i[np.where(y0<=y[i])]
         if not yy is None: i = i[np.where(y[i]<=yy)]
         
         if x0 is None: _x0 = np.min(np.append(x[i],_x0))
         if xx is None: _xx = np.max(np.append(x[i],_xx))
         if y0 is None: _y0 = np.min(np.append(y[i],_y0))
         if yy is None: _yy = np.max(np.append(y[i],_yy))

      if x0 is None and xx is None:
         x0, xx = _x0 - .03*(_xx-_x0), _xx + .03*(_xx-_x0)
      elif x0 is None: x0 = _x0 - .03*(xx-_x0)
      elif xx is None: xx = _xx + .03*(_xx-x0)

      if y0 is None and yy is None:
         y0, yy = _y0 - .05*(_yy-_y0), _yy + .05*(_yy-_y0)
      elif y0 is None: y0 = _y0 - .05*(yy-_y0)
      elif yy is None: yy = _yy + .05*(_yy-y0)
      
      if np.isfinite(x0+xx) and x0 < xx:
         ax.set_xlim(x0,xx)
      if np.isfinite(y0+yy) and y0 < yy:
         if self.currentfigure in ['Mag LC', 'contour', 'Summarize']: # invert y axis 
            ax.set_ylim(yy,y0)
         else:
            ax.set_ylim(y0,yy)         
      self.canvas.draw()
      
   def motion(self, e):
      x, y = e.xdata, e.ydata
      if x is not None or y is not None:
         self.Message('Position: {}, {}'.format(x, y))

   def move_from(self, e):
      ''' Remember previous coordinates for scrolling with the mouse '''
      ax = e.inaxes
      if not ax: return     
      # first click location
      self.movefromx, self.movefromy = e.xdata, e.ydata
      # first click time
      if not 'lastclicktime' in self.__dict__:
         self.lastclicktime = time.time()
         dt = 1
      else:
         dt = time.time() - self.lastclicktime
      if dt < .3: self.doubleclick = True         
      else:  self.doubleclick = False
      self.lastclicktime = time.time()
      # right/left click
      self.rightleftclick = e.button   
      
   def move_to(self, e):
      ''' Drag (move) canvas to the new position '''
      ax = e.inaxes
      if not ax: return
      
      if self.markerVar.get() == 1: self.move_to_marker(e)
      elif self.zoomVar.get() == 1: self.move_to_zoom(e)
      else:  self.move_to_pan(e)
      
   def move_to_marker(self, e, markersize=12):
      '''
      left:
         double click: pentagon
         single click: points
         holdon move: line
      right:
         double click: square
         single click: triangle
         holdon move: rectangle
      '''
      ax = e.inaxes
      if not ax: return
      if not 'movefromx' in self.__dict__: return
      # location shift
      dloc = np.sqrt( (e.xdata - self.movefromx)**2 + (e.ydata - self.movefromy) )**2      
      if self.rightleftclick == 1:  # left click
         if self.doubleclick: # double clock
            # pentagon
            self.ax.plot(e.xdata, e.ydata, color='r',
                         marker='p', fillstyle='none', markersize=markersize)
         else: # single click
            if dloc == 0:
               # point
               self.ax.plot(e.xdata, e.ydata, color='r', marker='o',
                            fillstyle='none', markersize=markersize)
            else:
               # line
               self.ax.plot([self.movefromx, e.xdata], [self.movefromy, e.ydata],
                            color='r', ls='--', fillstyle='none')         
      else:  # right click
         if self.doubleclick: # double clock
            self.ax.plot(e.xdata, e.ydata, color='r',
                     marker='s', fillstyle='none', markersize=markersize)
         else: # single click
            if dloc == 0:
               self.ax.plot(e.xdata, e.ydata, color='r',
                     marker='^', fillstyle='none', markersize=markersize)               
            else:
               # plot rectangle               
               self.ax.plot([self.movefromx, self.movefromx], [self.movefromy, e.ydata],
                            color='r', ls='-')
               self.ax.plot([self.movefromx, e.xdata], [self.movefromy, self.movefromy],
                            color='r', ls='-')
               self.ax.plot([self.movefromx, e.xdata], [e.ydata, e.ydata],
                            color='r', ls='-')
               self.ax.plot([e.xdata, e.xdata], [e.ydata, self.movefromy],
                            color='r', ls='-') 
      self.canvas.draw()
      
   def move_to_pan(self, e):
      '''
      left:
         double click: restart
         single click: -
         holdon move: pan
      right:
         double click: savefig
         single click: savefig
         holdon move:  pan
      '''
      ax = e.inaxes
      if not ax: return
      if not 'movefromx' in self.__dict__: return
      # location shift
      dloc = np.sqrt( (e.xdata - self.movefromx)**2 + (e.ydata - self.movefromy) )**2      
      if self.rightleftclick == 1:  # left click
         if self.doubleclick: # restart fig
            if self.currentfiguresize is not None:
               xlim = self.currentfiguresize[0]
               ylim = self.currentfiguresize[1]
               ax.set_xlim( xlim )
               ax.set_ylim( ylim )
         else: # single click
            if dloc == 0:
               return             
            else:
               # pan
               xs, ys = e.xdata - self.movefromx, e.ydata - self.movefromy
               x0,xx = ax.get_xlim()
               y0,yy = ax.get_ylim()
               self.set_limits(ax, x0=x0-xs, xx=xx-xs, y0=y0-ys, yy=yy-ys)         
      else:  # right click
         if self.doubleclick: # restart fig
            self.save_figure()
         else: # single click
            if dloc == 0:
               self.save_figure()
            else:
               # pan
               xs, ys = e.xdata - self.movefromx, e.ydata - self.movefromy
               x0,xx = ax.get_xlim()
               y0,yy = ax.get_ylim()
               self.set_limits(ax, x0=x0-xs, xx=xx-xs, y0=y0-ys, yy=yy-ys) 
      self.canvas.draw()

   def move_to_zoom(self, e):
      '''
      left:
         double click: restart
         single click: -
         holdon move: zoom
      right:
         double click: savefig
         single click: savefig
         holdon move: zoom
      '''
      ax = e.inaxes
      if not ax: return
      if not 'movefromx' in self.__dict__: return
      # location shift
      dloc = np.sqrt( (e.xdata - self.movefromx)**2 + (e.ydata - self.movefromy) )**2      
      if self.rightleftclick == 1:  # left click
         if self.doubleclick: # restart fig
            if self.currentfiguresize is not None:
               xlim = self.currentfiguresize[0]
               ylim = self.currentfiguresize[1]
               ax.set_xlim( xlim )
               ax.set_ylim( ylim )
         else: # single click
            if dloc == 0:
               return             
            else:
               # zoom
               if e.xdata < self.movefromx: ax.set_xlim(e.xdata,self.movefromx)
               else:  ax.set_xlim(self.movefromx,e.xdata)

               ymin, ymax = min(e.ydata, self.movefromy), max(e.ydata, self.movefromy)
               if self.currentfigure in ['Mag LC', 'contour', 'Summarize']: # invert y axis
                  ax.set_ylim(ymax, ymin)
               else:
                  ax.set_ylim(ymin, ymax)       
      else:  # right click
         if self.doubleclick: # restart fig
            self.save_figure()
         else: # single click
            if dloc == 0:
               self.save_figure()
            else:
               # zoom
               if e.xdata < self.movefromx: ax.set_xlim(e.xdata,self.movefromx)
               else:  ax.set_xlim(self.movefromx,e.xdata)

               ymin, ymax = min(e.ydata, self.movefromy), max(e.ydata, self.movefromy)
               if self.currentfigure in ['Mag LC', 'contour']: # invert y axis
                  ax.set_ylim(ymax, ymin)
               else:
                  ax.set_ylim(ymin, ymax)
      self.canvas.draw()
      
   def zoom(self, e):
      
      ax = e.inaxes
      if not ax: return            
      if not self.zoomable: return
      # x
      x,(x0,xx) = e.xdata, ax.get_xlim()
      if self.zoomable == 'right': x = x0
      
      if e.button == 2 and self.zoomable == 'horizontal': # center
         x0,xx = x-(xx-x0)/2, x+(xx-x0)/2
      elif e.button == 'down':
         x0,xx = x-0.8*abs(x-x0), x+0.8*abs(xx-x)
      elif e.button == 'up':
         x0,xx = x-1.2*abs(x-x0), x+1.2*abs(xx-x)

      # y
      y,(y0,yy) = e.ydata, ax.get_ylim()
      if self.zoomable == 'up': y = y0
      
      if e.button == 2 and self.zoomable == 'horizontal': # center
         y0,yy = y-(yy-y0)/2, y+(yy-y0)/2
      elif e.button == 'down':
         y0,yy = y-0.8*abs(y-y0), y+0.8*abs(yy-y)
      elif e.button == 'up':
         y0,yy = y-1.2*abs(y-y0), y+1.2*abs(yy-y)
         
      self.set_limits(ax, x0=x0, xx=xx, y0=y0, yy=yy)

   def add_subplot(self, position, zoomable=False, cursor=False):      
      if type(position) is SubplotSpec:         
         self.ax = self.fig.add_subplot(position)
      else:
         if type(position) is int: position = [position]         
         self.ax = self.fig.add_subplot(*position)         
      self.zoomable = zoomable
      if cursor:
         horizOn = True if cursor in ('horisontal','both') else False
         vertOn = True if cursor in ('vertical','both') else False
         cursor = Cursor(self.ax, horizOn=horizOn, vertOn=vertOn, useblit=True, color='black', lw=1)        
         self.cursor = cursor
      return self.ax
   
   def save_figure(self):
      
      filetypes = [
         ('Portable Network Graphics (*.png)', '*.png'),
         ('Encapsulated Postscript (*.eps)', '*.eps'),
         ('Portable Document Format', '*.pdf'),
         ('PGF code for LaTeX (*.pgf)', '*.pgf'),
         ('Postscript (*.ps)', '*.ps'),
         ('Raw RGBA bitmap (*.raw, *.rgba)', '*.raw, *.rgba'),
         ('Scalable Vector Graphics (*.svg, *.svgz)', '*.svg, *.vgz'),
      ]

      filename = filedialog.asksaveasfilename(parent=self, initialdir=os.getcwd(),
                                 title='Save the figure', filetypes=filetypes)

      if filename: self.fig.savefig(filename)       

   ######################################################   
   ######################################################                    
   def loadall(self, debug=True):
      objlist = []
      failed_objs = []
      with get_progress_bar(True, len(self.ztfp.meta.index)) as pbar:            
         for i, objid in enumerate(self.ztfp.meta.index):
            # load data
            if debug:
               loaded = self.ztfp.load_data(objid)
            else:
               try:
                  loaded = self.ztfp.load_data(objid)
               except:
                  failed_objs.append(objid)
                  continue
            if loaded:  objlist.append(objid)
            
            # update pbar
            pbar.update(1)
      
      if len(objlist) == 0:
         self.Message('No cached objects found')
         return
      self.Message('%i objects have been reloaded\n%i objects failed: %s'%
                   (len(objlist), len(failed_objs), failed_objs))
      
      # set objid
      self.setobjid(sorted(np.unique(objlist))[0])
      
      # update option menu
      menu = self.objlistMenu["menu"]
      menu.delete(0, "end")
      for t in sorted(np.unique(objlist)):
         menu.add_command(label=t,command=lambda k=t: self.setobjid(k))
         
      # renew header infos
      self.metainfo.set('Obj: %s \n%i of %i\n%i in total' %
                        (self.obj.get(), 0, len(objlist), len(self.ztfp.meta)))
      
   def saveall(self):
      
      #if self.clobberVar.get()==0: clobber = False
      #else: clobber = True
      clobber = True
      
      objid = self.obj.get()
      if objid in ['','None']:         
         objlist = list(self.ztfp.data.keys())
         
         for objid in objlist:
            saved = self.ztfp.save_data(objid, clobber=clobber)               
      else:
         saved = self.ztfp.save_data(objid, clobber=clobber)  
         
   def addobj(self):      
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: set object first')
         return
      objlist = list(self.ztfp.data.keys())
      if objid in objlist: # and self.clobberVar.get()==0:
         # skip
         self.Message('Warning: %s already added'%objid)
         return 
      
      self.ztfp.parse_meta_all(self.ztfp.kwargs, objid)
      
      # load data
      self.ztfp.load_data(objid)
      
      # if not saved, reload
      if not objid in self.ztfp.data or self.ztfp.kwargs['clobber']:
         par = dict()
         if objid in self.ztfp.params: par = self.ztfp.params[objid]
         
         # for each object                         
         self.ztfp.data[objid] = snerun.snobject(objid, z=self.ztfp.z, ra=self.ztfp.ra, dec=self.ztfp.dec,
                            mkwebv=self.ztfp.mkwebv, hostebv=self.ztfp.hostebv, sntype=self.ztfp.sntype,
                            dm=self.ztfp.dm, jdpeak=self.ztfp.jdpeak, **par)
         self.Message('%s was added'% (objid))         

      # append obj
      objlist.append(objid)
      
      # update option menu
      menu = self.objlistMenu["menu"]
      menu.delete(0, "end")
      for t in sorted(np.unique(objlist)):
         menu.add_command(label=t,command=lambda k=t: self.setobjid(k))
         
      # renew header infos      
      self.setobjid(objid)
      
   def delobj(self):
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: No objects left')
         return
      del self.ztfp.data[objid]
      
      objlist = list(self.ztfp.data.keys())
      self.Message('deleted %s, %i objects left '%(objid, len(objlist)))
      if len(objlist)>0: self.setobjid(objlist[0])
      else: self.setobjid(None)
      
      # update option menu
      menu = self.objlistMenu["menu"]
      menu.delete(0, "end")
      for t in sorted(np.unique(objlist)):
         menu.add_command(label=t,command=lambda k=t: self.setobjid(k)) 
      
   def setobjid(self, obj):
      self.obj.set(obj)
      # renew header info
      objlist = np.array(sorted(self.ztfp.data.keys()))
      if len(objlist) == 0:
         self.Message('Add some objects first')
         return
      if obj is not None:
         try:
            __ = np.where(objlist == obj)[0][0]
         except:
            return
      else:
         __ = 0
      self.metainfo.set('Obj: %s \n%i of %i\n%i in total' %
                        (self.obj.get(), __, len(self.ztfp.data), len(self.ztfp.meta)))
      
   def SearchObj(self):
      _str = self.searchFrameTxt.get()
      if _str == '':  # empty for tutorial
         self.Message('input %s' % searchnames)
         return
      objid = self.matchstr(_str, True)
      if objid is not None: self.setobjid(objid)
         
   def matchstr(self, _str, verbose=False):
      if len(self.ztfp.meta) == 0:
         msg = 'Warning:\n\nno meta, parse it or manually input parameters'
         flag = None
      else:
         candlist = []
         with get_progress_bar(verbose, len(self.ztfp.meta.index)) as pbar:             
            for i, objid in enumerate(self.ztfp.meta.index):
               if  _str in objid:
                  candlist.append(objid)
                  continue
               for searchname in searchnames:                  
                  if searchname in self.ztfp.meta.keys():                   
                     kid = str(self.ztfp.meta[searchname][i])
                     if _str in kid: candlist.append(objid)                        
               pbar.update(1)
         
         if len(candlist) == 0:
            msg = 'Warning:\n\nno objects matched to %s'%_str
            flag = None
         elif len(candlist) == 1:
            msg = 'Nice:\n\n%s was matched successfully as %s\n\nInfo:\n\n%s'%\
               (_str,candlist[0], tabulate(self.ztfp.meta.query('%s=="%s"'%(snelist_pars['idkey'], candlist[0])), headers = 'keys', tablefmt = 'grid'))
            flag = candlist[0]
         elif len(candlist) > 1:
            msg = 'Warning:\n\nMore than one object was matched with %s, specify one from below:\n\n%s'%\
               (_str,candlist)
            flag = None        
      if verbose:  self.Message(msg, False)
      return flag                          
      
   def runtasks(self):
      if self.optvar.get() is None:
         self.Message('No task specified')
         return
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: set object first')
         return      
      objlist = list(self.ztfp.data.keys())
      if objid not in objlist:
         self.Message('Warning: add object first')
         return      
      if self.optvar.get() == 'Data preparation':
         if not 'datatop' in self.__dict__: self.data_popup(objid)
         elif self.datatop.winfo_exists() == 0: self.data_popup(objid)
         else: self.Message('Parse photometric/spectral topup window already exists')
      if self.optvar.get() == 'Multiband LC fits':
         if not 'lc' in self.ztfp.data[objid].__dict__:
            self.Message('Parse LC data for %s first before fitting'%objid)
            return
         if not 'datatop1' in self.__dict__: self.data_popup1(objid)
         elif self.datatop1.winfo_exists() == 0: self.data_popup1(objid)
         else: self.Message('Fit photometric/spectral topup window already exists')
      if self.optvar.get() == 'SED fits':
         if not 'lc' in self.ztfp.data[objid].__dict__:
            self.Message('Parse LC data for %s first before fitting'%objid)
            return
         if not 'bbtop' in self.__dict__: self.bb_popup(objid)
         elif self.bbtop.winfo_exists() == 0: self.bb_popup(objid)
         else: self.Message('Parse Colours/Bolometric topup window already exists')
      if self.optvar.get() == 'Bolometric LC fits':         
         if not 'bbtop1' in self.__dict__: self.bb_popup1(objid)
         elif self.bbtop1.winfo_exists() == 0: self.bb_popup1(objid)
         else: self.Message('Fit Bolometric topup window already exists')          
      if self.optvar.get() == 'Spectral line fits':
         if not 'spectop' in self.__dict__: self.spectrafit_popup(objid)
         elif self.spectop.winfo_exists() == 0: self.spectrafit_popup(objid)
         else: self.Message('Spectral Line Fit topup window already exists')           
         
      self.Message('task %s complete' % self.optvar.get())
      # reset action var
      #self.optvar.set(None)
      
   def runplots(self):
      if self.plotvar.get() is None:
         self.Message('Set proper figstyle first')
         return     
      objid = self.obj.get()
      objlist = list(self.ztfp.data.keys())
      if objid == 'None':
         self.Message('Warning: set object first')
         return
      if objid not in objlist:
         self.Message('Warning: add %s first'%objid)
         return
      # reset canvas
      try:    self.fig.clear(False)
      except: pass
      # change current figrue
      self.currentfigure = self.plotvar.get()
      # plot
      if self.plotvar.get() == 'Summarize Plot':         
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].summarize_plot(ax=ax, savefig=True)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show summary lcs for %s'%objid)
         
      if self.plotvar.get() == 'Flux LC':
         if not 'lc' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse lc for %s first'%objid)
            return
         ax = self.ztfp.data[objid].ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax(show_title=False, ax_xstyle='jd', ax_ystyle='original')
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show Flux lcs for %s'%objid)

         if self.fitVar.get() == 1:
            if not 'fluxtop' in self.__dict__: self.fluxfit_popup(objid)
            elif self.fluxtop.winfo_exists() == 0: self.fluxfit_popup(objid)
            else: self.Message('Flux plot window already exists')
         
      if self.plotvar.get() == 'Mag LC':
         if not 'lc' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse lc for %s first'%objid)
            return          
         ax = self.ztfp.data[objid].ax2 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax2(show_legend=True, show_texp=True)                                    
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show Mag lcs for %s'%objid)

      if self.plotvar.get() == 'colour': 
         ax = self.ztfp.data[objid].ax3 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax3(show_legend=True, show_texp=True, show_fits=True, show_gp=True, ylabel_2right=False)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show colour curve for %s'%objid)          

      if self.plotvar.get() == 'luminosity':
         if not 'mbol' in self.ztfp.data[objid].__dict__.keys() and \
            not 'mbolbb' in self.ztfp.data[objid].__dict__.keys():         
            self.Message('Warning: parse bolometric LCs for %s first'%objid)
            return
         ax = self.ztfp.data[objid].ax4 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         make_bol = []
         if 'mbol' in self.ztfp.data[objid].__dict__: make_bol.append('lyman')
         if 'mbolbb' in self.ztfp.data[objid].__dict__: make_bol.append('bb')
         if 'mbolspec' in self.ztfp.data[objid].__dict__: make_bol.append('spec')
         if len(make_bol) == 0:
            self:message('no bolometric LCs found, try get some first')
            return
         self.ztfp.data[objid]._ax4(
            make_bol=make_bol, show_legend=True, show_texp=True, show_fits=True
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show bolometric LCs for %s'%objid)

         if self.fitVar.get() == 1:
            if not 'mfittop' in self.__dict__: self.mbolfit_popup()
            elif self.mfittop.winfo_exists() == 0: self.mbolfit_popup()
            else: self.Message('model mbol generator window already exists')
         
      if self.plotvar.get() == 'LC baseline':
         if not 'lc' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse lc for %s first'%objid)
            return         
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         baseline = self.ztfp.data[objid].calibrate_baseline(
            ax=ax, key='fcqfid', source='ztffp',
            xmin=-100, xmax=-20, ax_xlim=None, ax_ylim=None
         )
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('calibrate baseline for %s %s' % (objid,baseline))

         if self.fitVar.get() == 1:
            if not 'baselinetop' in self.__dict__: self.baselinetop_popup(objid)
            elif self.baselinetop.winfo_exists() == 0: self.baselinetop_popup(objid)
            else: self.Message('Baseline popup window already exists')
            
      if self.plotvar.get() == 'spectra':
         if not 'spec' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse spectra for %s first'%objid)
            return
         ax = self.ztfp.data[objid].ax1 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax1(text_type='phase')
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show spectra for %s'%objid)         
         
      if self.plotvar.get() == 'Finder':
         if not 'ra' in self.ztfp.data[objid].__dict__.keys() or \
            not 'dec' in self.ztfp.data[objid].__dict__.keys():         
            self.Message('Warning: parse coordinate for %s first'%objid)
            return
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         fcls = image_tool.handle_image(self.ztfp.data[objid].ra, self.ztfp.data[objid].dec,
                           img='%s/images/%s.fits' % (LOCALSOURCE, objid ), ax=ax, fov=100,
                           source='ps1', filt='r', clobber=False)
         fcls.run()
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show stamp for %s'%objid)
         
      if self.plotvar.get() == 'SNe scatter':
         if not 'disttop' in self.__dict__: self.dist_popup()
         elif self.disttop.winfo_exists() == 0: self.dist_popup()
         else: self.Message('Feature distribution topup window already exists')

      self.Message('Plot %s complete' % self.plotvar.get())
      # reset action var
      #self.plotvar.set(None)
      
   def checkdata(self, objid=None, checkitems=None):
      if objid is None:  objid = self.obj.get()
      objlist = list(self.ztfp.data.keys())
      if objid == 'None':  return 'Warning: set object first'         
      elif objid not in objlist:  return 'Warning: add %s first'%objid         
      elif checkitems is None:  return 'Warning: set checkitems'         
      msg = ''
      if 'data' in checkitems:
         if 'lc' in self.ztfp.data[objid].__dict__:
            lc = self.ztfp.data[objid].lc
            msg += 'LC: \n'
            for s in np.unique(lc['source']): msg += ' - %s\n' % (s)
         else: msg += 'LC: - \n'
         if 'spec' in self.ztfp.data[objid].__dict__:
            msg += 'spectra:\n  %i epochs \n' % len(self.ztfp.data[objid].spec.data)
         else:
            msg += 'spectra: - \n'
      if 'fits' in checkitems:
         if 'gpcls' in self.ztfp.data[objid].__dict__:
            msg += 'GP: Done \n'
         else:
            msg += 'GP: - \n'
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            msg += 'Fits: \n'
            for engine in self.ztfp.data[objid].fitcls:
               if engine not in ['multiband_main', 'multiband_early']: continue
               msg += ' - %s\n' % (engine)
               for model in self.ztfp.data[objid].fitcls[engine]:                  
                  msg += '    %s \n' % model      
         else:
            msg += 'Fits: - \n'         
      if 'bc' in checkitems:
         if 'colors' in self.ztfp.data[objid].__dict__:         
            msg += 'colors: Done \n'
         else:
            msg += 'colors: - \n'
         if 'mbol' in self.ztfp.data[objid].__dict__:
            msg += 'Lyman mbol: Done \n'
         else:
            msg += 'Lyman mbol: - \n'
         if 'mbolbb' in self.ztfp.data[objid].__dict__:
            msg += 'blackbody mbol: Done \n'
         else:
            msg += 'blackbody mbol: - \n'
         if 'mbolspec' in self.ztfp.data[objid].__dict__:
            msg += 'spectra mbolspec: Done \n'
         else:
            msg += 'spectra mbolspec: - \n'
      if 'bolfit' in checkitems:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            msg += 'Fits: \n'
            for engine in self.ztfp.data[objid].fitcls:
               if engine not in ['bol_main', 'bol_early', 'bol_tail', 'bol_full']: continue
               msg += ' - %s\n' % (engine)               
         else:
            msg += 'Fits: - \n'
      if 'sed' in checkitems:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            msg += 'Fits: \n'
            for engine in self.ztfp.data[objid].fitcls:
               if engine not in ['specline', 'specv_evolution']: continue               
               for source in self.ztfp.data[objid].fitcls[engine]:                  
                  msg += ' %s \n'%source
         else:
            msg += 'Fits: - \n'         
      if 'sedfits' in checkitems:
         if 'fitcls' in self.ztfp.data[objid].__dict__:            
            msg += 'SED Fits: \n'
            if 'sed' in self.ztfp.data[objid].fitcls:
               bb, spec = 0, 0
               for source in self.ztfp.data[objid].fitcls['sed']:
                  if 'bb' in source: bb += 1
                  if 'spec' in source: spec += 1
               if bb > 0: msg += ' bb: %i' % bb
               if spec > 0: msg += ' spec: %i' % spec
         else:
            msg += 'Fits: - \n'
      return msg

   ### popup for scatter ###
   def dist_popup(self):
      if self.clobberVar.get()==0: clobber = False
      else: clobber = True
      
      # pop up window
      self.disttop = Tk.Toplevel(self.parent)
      self.disttop.geometry("500x300+50+50")      
      self.disttop.title('Distribute features for %i objects'%len(self.ztfp.data))
      
      # subset frame
      labelFrame = Tk.LabelFrame(self.disttop, text="")
      labelFrame.grid(row=0,column=0,sticky=Tk.NW,rowspan=1)
      self.subsetinfo = Tk.StringVar()
      self.subsetinfo.set('Data subset:\n%s' % ( '\n'.join(self.ztfp.dset.keys())))
      tb = Tk.Label(labelFrame,textvariable=self.subsetinfo,fg="red")      
      tb.grid(row=0,column=0)
      
      # parameter frame
      labelFrame = Tk.LabelFrame(self.disttop, text="")
      labelFrame.grid(row=0,column=2,sticky=Tk.NW,rowspan=1)
      self.subsetinfo1 = Tk.StringVar()
      par = []
      for _par in self.ztfp.pars:
         par.append(self.ztfp.get_par(None, returnname=True, **_par))         
      self.subsetinfo1.set('Parameters:\n%s' % ( '\n'.join(par)))
      tb = Tk.Label(labelFrame,textvariable=self.subsetinfo1,fg="red")
      tb.grid(row=0,column=0)
      
      def subsetpop():
         if not 'subsettop' in self.__dict__: self.subset_popup()
         elif self.subsettop.winfo_exists() == 0: self.subset_popup()
         else: self.Message('Add subset topup window already exists')
      data1button = Tk.Button(self.disttop,text="Popup",command=subsetpop)
      data1button.grid(row=0, column=1)
      
      def scatterppop():
         if not 'scatterptop' in self.__dict__: self.scatterp_popup()
         elif self.scatterptop.winfo_exists() == 0: self.scatterp_popup()
         else: self.Message('Add parameter topup window already exists')
      data1button = Tk.Button(self.disttop,text="Popup",command=scatterppop)
      data1button.grid(row=0, column=3)
      
      # 1d            
      Frame = Tk.LabelFrame(self.disttop, text="1D histogram")
      Frame.grid(row=1,column=0,sticky=Tk.W)
      
      l = Tk.Label(Frame,text='x=',fg="black")      
      l.grid(row=0,column=0)
      
      self.pval = Tk.StringVar()
      self.pval.set( None )
      self.scatter1dMenu = Tk.OptionMenu(Frame, self.pval, *[None])
      self.scatter1dMenu.config(justify=Tk.LEFT,width=8)
      self.scatter1dMenu.grid(row=0,column=1, sticky=Tk.W)

      l = Tk.Label(Frame,text='style=',fg="black")      
      l.grid(row=1,column=0)
      
      self.pstval = Tk.StringVar()
      self.pstval.set( 'pdf' )
      Menu = Tk.OptionMenu(Frame, self.pstval, *['pdf', 'cdf'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=1,column=1, sticky=Tk.W)
      
      # 2d
      Frame = Tk.LabelFrame(self.disttop, text="2D scatter")
      Frame.grid(row=2,column=0,sticky=Tk.W)
      
      l = Tk.Label(Frame,text='x=',fg="black")      
      l.grid(row=0,column=0)
      
      self.pval1 = Tk.StringVar()
      self.pval1.set( None )
      self.scatter2dxMenu = Tk.OptionMenu(Frame, self.pval1, *[None])
      self.scatter2dxMenu.config(justify=Tk.LEFT,width=8)
      self.scatter2dxMenu.grid(row=0,column=1, sticky=Tk.W)
      
      l = Tk.Label(Frame,text='y=',fg="black")      
      l.grid(row=1,column=0)
      
      self.pval2 = Tk.StringVar()
      self.pval2.set( None )
      self.scatter2dyMenu = Tk.OptionMenu(Frame, self.pval2, *[None])
      self.scatter2dyMenu.config(justify=Tk.LEFT,width=8)
      self.scatter2dyMenu.grid(row=1,column=1, sticky=Tk.W)  

      # nd            
      Frame = Tk.LabelFrame(self.disttop, text="nD Contours")
      Frame.grid(row=3,column=0,sticky=Tk.W)
      
      l = Tk.Label(Frame,text='dataset=',fg="black")      
      l.grid(row=0,column=0)
      
      self.pdataset = Tk.StringVar()
      self.pdataset.set( None )
      self.scatterndMenu = Tk.OptionMenu(Frame, self.pdataset, *[None])
      self.scatterndMenu.config(justify=Tk.LEFT,width=8)
      self.scatterndMenu.grid(row=0,column=1, sticky=Tk.W)
      
      # show all
      Frame = Tk.LabelFrame(self.disttop, text="show all")
      Frame.grid(row=4,column=0,columnspan=3,sticky=Tk.W)
      
      l = Tk.Label(Frame,text='figure=',fg="black")      
      l.grid(row=0,column=0)
      
      self.scatterfigure = Tk.StringVar()
      self.scatterfigure.set( 'flux' )
      Menu = Tk.OptionMenu(Frame, self.scatterfigure,
         *['flux', 'mag', 'colour', 'luminosity', 'velocity'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=1, sticky=Tk.W)

      l = Tk.Label(Frame,text='dataset=',fg="black")      
      l.grid(row=0,column=2)
      
      self.showalldataset = Tk.StringVar()
      self.showalldataset.set( None )
      self.showallMenu = Tk.OptionMenu(Frame, self.showalldataset, *[None])
      self.showallMenu.config(justify=Tk.LEFT,width=8)
      self.showallMenu.grid(row=0,column=3, sticky=Tk.W)
      
      l = Tk.Label(Frame,text='show data=',fg="black")      
      l.grid(row=1,column=0)
      
      self.showdatapoint = Tk.StringVar()
      self.showdatapoint.set( 'True' )
      Menu = Tk.OptionMenu(Frame, self.showdatapoint, *['True', 'False'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=1,column=1, sticky=Tk.W)  
      
      l = Tk.Label(Frame,text='show fit=',fg="black")      
      l.grid(row=1,column=2)
      
      self.showdatafit = Tk.StringVar()
      self.showdatafit.set( 'False' )
      Menu = Tk.OptionMenu(Frame, self.showdatafit, *['True', 'False'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=1,column=3, sticky=Tk.W)

      l = Tk.Label(Frame,text='show gp=',fg="black")      
      l.grid(row=2,column=0)
      
      self.showdatagp = Tk.StringVar()
      self.showdatagp.set( 'False' )
      Menu = Tk.OptionMenu(Frame, self.showdatagp, *['True', 'False'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=2,column=1, sticky=Tk.W)

      l = Tk.Label(Frame,text='fit error=',fg="black")      
      l.grid(row=2,column=2)
      
      self.showdatafiterr = Tk.StringVar()
      self.showdatafiterr.set( 'False' )
      Menu = Tk.OptionMenu(Frame, self.showdatafiterr, *['True', 'False'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=2,column=3, sticky=Tk.W)
      
      l = Tk.Label(Frame,text='filters=',fg="black")      
      l.grid(row=3,column=0)
      
      self.showfilterentry = Tk.Entry(Frame, justify=Tk.LEFT, width=8)      
      self.showfilterentry.insert(0, 'r')
      self.showfilterentry.grid(row=3, column=1, sticky=Tk.W)      
      
      # table      
      Frame = Tk.LabelFrame(self.disttop, text="Make table")
      Frame.grid(row=1,column=2,rowspan=3,sticky=Tk.W)
      
      l = Tk.Label(Frame,text='parameter=',fg="black")      
      l.grid(row=0,column=0)
      
      self.tablepentry = Tk.Entry(Frame, justify=Tk.LEFT, width=8)      
      self.tablepentry.insert(0, [])
      self.tablepentry.grid(row=0, column=1, sticky=Tk.W)

      l = Tk.Label(Frame,text='dataset=',fg="black")      
      l.grid(row=1,column=0)
      
      self.tabdataset = Tk.StringVar()
      self.tabdataset.set( None )
      self.stabMenu = Tk.OptionMenu(Frame, self.pdataset, *[None])
      self.stabMenu.config(justify=Tk.LEFT,width=8)
      self.stabMenu.grid(row=1,column=1, sticky=Tk.W)
      
      l = Tk.Label(Frame,text='style=',fg="black")      
      l.grid(row=2,column=0)
      
      self.tabstyle = Tk.StringVar()
      self.tabstyle.set( 'latex' )
      Menu = Tk.OptionMenu(Frame, self.tabstyle, *['latex', 'space', 'comma'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=2,column=1,sticky=Tk.W)
      
      # plots
      def run1d():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         p = self.pval.get()
         if p is None or p == 'None':
            self.Message('add a parameter first')
            return         
         self.currentfigure = 'scatter'
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         index = None
         for n, par in enumerate(self.ztfp.pars):
            parname = self.ztfp.get_par(None, returnname=True, **par)
            if parname == p: index = n
         if index is None: return
         self.ztfp.show1d(index=index, style=self.pstval.get())
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, parname)
         self.set_limits( ax )
         self.Message('1d distribution plots')         
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=run1d)
      inputFrameBut.grid(row=1, column=1, sticky=Tk.SW)
      
      def run2d():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         p1 = self.pval1.get()
         p2 = self.pval2.get()
         if (p1 is None or p1 == 'None') or (p2 is None or p2 == 'None'):
            self.Message('add a parameter first')
            return
         self.currentfigure = 'scatter'
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         index1, index2 = None, None
         for n, par in enumerate(self.ztfp.pars):
            parname = self.ztfp.get_par(None, returnname=True, **par)
            if parname == p1: index1, p1n = n, parname
            if parname == p2: index2, p2n = n, parname
         if index1 is None or index2 is None: return         
         self.ztfp.show2d(index1,index2)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, '%s vs %s'%(p1n,p2n))
         self.set_limits( ax )
         self.Message('2d distribution plots')
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=run2d)
      inputFrameBut.grid(row=2, column=1, sticky=Tk.SW)
      
      def runnd():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         syntax = self.pdataset.get()
         if syntax is None or syntax == 'None':
            self.Message('add a dataset first')
            return
         self.currentfigure = 'scatter'
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')        
         self.ztfp.shownd(syntax=syntax)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, 'scatter for dataset: %s'%syntax)
         self.set_limits( ax )
         self.Message('nd corner plots')
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=runnd)
      inputFrameBut.grid(row=3, column=1, sticky=Tk.SW)

      def runall():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         fig = self.scatterfigure.get()
         syntax = self.showalldataset.get()
         if syntax is None or syntax == 'None':
            self.Message('add a dataset first')
            return
         show_data = eval(self.showdatapoint.get())
         show_fits = eval(self.showdatafit.get())
         show_gp = eval(self.showdatagp.get())
         plot_bands = self.showfilterentry.get().split()
         show_fit_error = eval(self.showdatafiterr.get())
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         if fig == 'mag':
            self.currentfigure = 'Mag LC'
            self.ztfp.showax2(               
               syntax=syntax, plot_bands=plot_bands, show_data=show_data,
               show_fits=show_fits, show_gp=show_gp, show_fit_error=show_fit_error,
            ) 
         elif fig == 'flux':
            self.currentfigure = 'Flux LC'
            self.ztfp.showax(
               syntax=syntax, plot_bands=plot_bands, show_data=show_data,
               show_fits=show_fits, show_gp=show_gp, show_fit_error=show_fit_error
            )  
         elif fig == 'colour':
            self.currentfigure = 'colour curve'
            self.ztfp.showax3(
               syntax=syntax, show_data=show_data,
               show_fits=show_fits, show_gp=show_gp, show_fit_error=show_fit_error
            ) 
         elif fig == 'luminosity':
            self.currentfigure = 'luminosity LC'
            self.ztfp.showax4(
               syntax=syntax, show_data=show_data,
               show_fits=show_fits, show_fit_error=show_fit_error
            ) 
         elif fig == 'velocity':
            self.currentfigure = 'velocity curve'
            self.ztfp.showax6(
               syntax=syntax, show_data=show_data,
               show_fits=show_fits, show_fit_error=show_fit_error
            ) 
         else:
            self.Message('Error: wrong figure style')
            return
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, '%s'%self.currentfigure)
         self.set_limits( ax )
         self.Message('show %s plots for all objects'%self.currentfigure)
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=runall)
      inputFrameBut.grid(row=4, column=3, sticky=Tk.W) 
      
      def maketable():
         tablepars = self.tablepentry.get().split()
         style = self.tabstyle.get()
         syntax = self.tabdataset.get()
         table = self.ztfp.table(
            syntax=syntax, tablepars=tablepars, tablename=None, style=style,
         )
         self.Message(table, disable=False)
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=maketable)
      inputFrameBut.grid(row=3, column=3, sticky=Tk.W)      
      
   def scatterp_popup(self):
      # pop up window
      self.scatterptop = Tk.Toplevel(self.parent)
      self.scatterptop.geometry("300x400+50+50")      
      self.scatterptop.title('Add parameter')
      
      # pars
      Frame = self.scatterptop
      filters = list(central_wavelengths.keys())
      
      rownumber = 0
      l = Tk.Label(Frame,text='parname:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterparname = Tk.StringVar()
      self.scatterparname.set( snerun.snelist.parlist[0] )      
      Menu = Tk.OptionMenu(Frame, self.scatterparname, *snerun.snelist.parlist)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)

      rownumber = 1      
      l = Tk.Label(Frame,text='filter 1:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterfilter1 = Tk.StringVar()
      self.scatterfilter1.set( 'g' )      
      Menu = Tk.OptionMenu(Frame, self.scatterfilter1, *filters)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)
      
      rownumber = 2
      l = Tk.Label(Frame,text='filter 2:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterfilter2 = Tk.StringVar()
      self.scatterfilter2.set( 'r' )      
      Menu = Tk.OptionMenu(Frame, self.scatterfilter2, *filters)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)

      rownumber = 3
      l = Tk.Label(Frame,text='peak with:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterpeakwith = Tk.StringVar()
      self.scatterpeakwith.set( 'gp' )      
      Menu = Tk.OptionMenu(Frame, self.scatterpeakwith, *['gp', 'multiband', 'bol'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)

      rownumber = 4
      l = Tk.Label(Frame,text='explosion with:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterexpwith = Tk.StringVar()
      self.scatterexpwith.set( 'mid' )      
      Menu = Tk.OptionMenu(Frame, self.scatterexpwith, *['mid', 'pl', 'bol'])
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)
      
      rownumber = 5
      l = Tk.Label(Frame,text='phase 1:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterphase1 = Tk.Entry(Frame, justify=Tk.LEFT) 
      self.scatterphase1.insert(0, 0)
      self.scatterphase1.config(justify=Tk.LEFT,width=12)
      self.scatterphase1.grid(row=rownumber, column=1)
      
      rownumber = 6
      l = Tk.Label(Frame,text='phase 2:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterphase2 = Tk.Entry(Frame, justify=Tk.LEFT)      
      self.scatterphase2.insert(0, 0)
      self.scatterphase2.config(justify=Tk.LEFT,width=12)
      self.scatterphase2.grid(row=rownumber, column=1)
      
      rownumber = 7
      l = Tk.Label(Frame,text='interpolation:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterinterpolation = Tk.StringVar()
      cinterp = [None, 'bin','gp','fit']
      self.scatterinterpolation.set( cinterp[0] )
      Menu = Tk.OptionMenu(Frame, self.scatterinterpolation, *cinterp)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)
      
      rownumber = 8
      l = Tk.Label(Frame,text='corr mkw:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scattercorrmkwebv = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(Frame, text="",variable=self.scattercorrmkwebv) 
      cButton.grid(row=rownumber,column=1,sticky=Tk.W)
      
      rownumber = 9
      l = Tk.Label(Frame,text='corr host:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scattercorrhostebv = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(Frame, text="",variable=self.scattercorrhostebv) 
      cButton.grid(row=rownumber,column=1,sticky=Tk.W)

      engines = [None] + list(get_engine().keys())
      
      rownumber = 11
      l = Tk.Label(Frame,text='model:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scattermodel = Tk.StringVar()
      opts = [None] + get_model()[engines[1]]      
      self.scattermodel.set( opts[0] )
      self.modelMenu = Tk.OptionMenu(Frame, self.scattermodel, *opts)
      self.modelMenu.config(justify=Tk.LEFT,width=8)
      self.modelMenu.grid(row=rownumber,column=1, sticky=Tk.W)
      
      rownumber = 10
      l = Tk.Label(Frame,text='engine:',fg="black")      
      l.grid(row=rownumber,column=0)
      self.scatterengine = Tk.StringVar()      
      self.scatterengine.set( engines[0] )
      def updatemodel1(k):
         self.scattermodel.set( k )
      def updatemodel(engine):
         opts = get_model()[engine]
         # update option menu
         self.scattermodel.set( get_model()[engine][0] )
         menu = self.modelMenu["menu"]
         menu.delete(0, "end")
         for t in opts:
            menu.add_command(label=t,command=lambda k=t: updatemodel1(k)) 
      Menu = Tk.OptionMenu(Frame, self.scatterengine, *engines,
            command=lambda k=self.scatterengine.get(): updatemodel(k))
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=rownumber,column=1, sticky=Tk.W)      
      
      def addpar():
         if self.scattercorrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.scattercorrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         if self.scatterengine.get() in ['', 'None']:
            scatterengine = None
         else:
            scatterengine = self.scatterengine.get()
         if self.scattermodel.get() in ['', 'None']:
            scattermodel = None
         else:
            scattermodel = self.scattermodel.get()
         if self.scatterinterpolation.get() in ['', 'None']:
            cinterp = None
         else:
            cinterp = self.scatterinterpolation.get()
         kwargs = {
            'parname'      : self.scatterparname.get(),
            'filt1'        : self.scatterfilter1.get(),
            'filt2'        : self.scatterfilter2.get(),
            'quant'        : snobject_pars['quantile'],
            'peakwith'     : self.scatterpeakwith.get(),
            'expwith'      : self.scatterexpwith.get(),
            'interpolation': cinterp,
            'phase1'       : float(self.scatterphase1.get()),
            'phase2'       : float(self.scatterphase2.get()),
            'corr_mkw'     : corr_mkw,
            'corr_host'    : corr_host,
            'index'        : 0,
            'engine'       : scatterengine,
            'model'        : scattermodel,
            'source'       : None,
         }
         self.ztfp.add_parameter(**kwargs)
         par = []
         for _par in self.ztfp.pars:
            par.append(self.ztfp.get_par(None, returnname=True, **_par))         
         self.subsetinfo1.set('Parameters:\n%s' % ( '\n'.join(par)))
         
         # update menu
         def updateval(tokey, val): tokey.set(val)

         # 1d
         self.pval.set( par[0] )
         menu = self.scatter1dMenu["menu"]
         menu.delete(0, "end")
         for t in par:
            menu.add_command(label=t,command=lambda k=t: updateval(self.pval, k))   
         # 2d x
         self.pval1.set( par[0] )
         menu = self.scatter2dxMenu["menu"]
         menu.delete(0, "end")
         for t in par:
            menu.add_command(label=t,command=lambda k=t: updateval(self.pval1, k))   
         # 2d y
         self.pval2.set( par[0] )
         menu = self.scatter2dyMenu["menu"]
         menu.delete(0, "end")
         for t in par:
            menu.add_command(label=t,command=lambda k=t: updateval(self.pval2, k))
         # table
         self.tablepentry.delete(0, "end")
         self.tablepentry.insert(0, par )         
      inputFrameBut = Tk.Button(Frame, text="Add", command=addpar)
      inputFrameBut.grid(row=12, column=0, sticky=Tk.W)
      
      def delpar():                  
         if len(self.ztfp.pars) == 0: return
         self.ztfp.pars.pop()         
         par = []
         for _par in self.ztfp.pars:
            par.append(self.ztfp.get_par(None, returnname=True, **_par))         
         self.subsetinfo1.set('Parameters:\n%s' % ( '\n'.join(par)))                  
         # update menu
         # 1d
         if len(self.ztfp.pars) == 0: self.pval.set( None ) 
         else: self.pval.set( par[0] )
         menu = self.scatter1dMenu["menu"]
         menu.delete(0, "end")
         if len(self.ztfp.pars) > 0: 
            for t in par:
               menu.add_command(label=t,command=lambda k=t: updateval(self.pval, k))   
         # 2d x
         if len(self.ztfp.pars) == 0: self.pval1.set( None ) 
         else: self.pval1.set( par[0] )
         menu = self.scatter2dxMenu["menu"]
         menu.delete(0, "end")
         if len(self.ztfp.pars) > 0: 
            for t in par:
               menu.add_command(label=t,command=lambda k=t: updateval(self.pval1, k))   
         # 2d y
         if len(self.ztfp.pars) == 0: self.pval2.set( None ) 
         else: self.pval2.set( par[0] )
         menu = self.scatter2dyMenu["menu"]
         menu.delete(0, "end")
         if len(self.ztfp.pars) > 0: 
            for t in par:
               menu.add_command(label=t,command=lambda k=t: updateval(self.pval2, k))
         # table
         self.tablepentry.delete(0, "end")
         self.tablepentry.insert(0, par ) 
         
      inputFrameBut = Tk.Button(Frame, text="Del", command=delpar)
      inputFrameBut.grid(row=12, column=1, sticky=Tk.W)
      
   def subset_popup(self):
      # pop up window
      self.subsettop = Tk.Toplevel(self.parent)
      self.subsettop.geometry("300x400+50+50")      
      self.subsettop.title('Add data subset')
      
      # syntax      
      Frame = Tk.LabelFrame(self.subsettop, text="add subset", width=6)
      Frame.grid(row=0,column=0,sticky=Tk.NW)
      
      l = Tk.Label(Frame,text='astype:',fg="black")      
      l.grid(row=0,column=0) 
      self.astypesyntax = Tk.Entry(Frame, justify=Tk.LEFT)
      syntax = 'z:float'
      self.astypesyntax.insert(0, syntax)
      self.astypesyntax.grid(row=0, column=1)
      
      l = Tk.Label(Frame,text='syntax:',fg="black")      
      l.grid(row=1,column=0)      
      self.searchsyntax = Tk.Entry(Frame, justify=Tk.LEFT)
      syntax = 'None'   
      self.searchsyntax.insert(0, syntax)
      self.searchsyntax.grid(row=1, column=1)
      
      # plotter
      Frame = Tk.LabelFrame(self.subsettop, text="subset plotter", width=6)
      Frame.grid(row=1,column=0,sticky=Tk.NW)
      
      l = Tk.Label(Frame,text='color:',fg="black")      
      l.grid(row=0,column=0)      
      self.subsetcolor = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetcolor.insert(0, '"k"')
      self.subsetcolor.grid(row=0, column=1)
      
      l = Tk.Label(Frame,text='nbin:',fg="black")      
      l.grid(row=1,column=0)
      self.subsetnbin = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetnbin.insert(0, '10')
      self.subsetnbin.grid(row=1, column=1)

      l = Tk.Label(Frame,text='fontsize:',fg="black")      
      l.grid(row=2,column=0)
      self.subsetfontsize = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetfontsize.insert(0, '12')
      self.subsetfontsize.grid(row=2, column=1)

      l = Tk.Label(Frame,text='label:',fg="black")      
      l.grid(row=3,column=0)
      self.subsetlabel = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetlabel.insert(0, 'None')
      self.subsetlabel.grid(row=3, column=1)

      l = Tk.Label(Frame,text='linestyle:',fg="black")      
      l.grid(row=4,column=0)
      self.subsetls = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetls.insert(0, '"-"')
      self.subsetls.grid(row=4, column=1)

      l = Tk.Label(Frame,text='marker:',fg="black")      
      l.grid(row=5,column=0)
      self.subsetmarker = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetmarker.insert(0, '"o"')
      self.subsetmarker.grid(row=5, column=1)
      
      l = Tk.Label(Frame,text='markersize:',fg="black")      
      l.grid(row=6,column=0)
      self.subsetmarkersize = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetmarkersize.insert(0, '12')
      self.subsetmarkersize.grid(row=6, column=1)

      l = Tk.Label(Frame,text='fillstyle:',fg="black")      
      l.grid(row=7,column=0)
      self.subsetfillstyle = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetfillstyle.insert(0, 'None')
      self.subsetfillstyle.grid(row=7, column=1)
      
      l = Tk.Label(Frame,text='labelpad:',fg="black")      
      l.grid(row=8,column=0)
      self.subsetlabelpad = Tk.Entry(Frame, justify=Tk.LEFT)   
      self.subsetlabelpad.insert(0, '12')
      self.subsetlabelpad.grid(row=8, column=1)
      
      def addsubset():
         syntax = self.searchsyntax.get()         
         if syntax in ['', 'None']: syntax = None
         astype = dict()
         for k in  self.astypesyntax.get().split(','):
            if ':' not in k: continue
            if len(k.split(':')) != 2: continue
            astype[k.split(':')[0]] = eval(k.split(':')[1])
         if len(astype) == 0: astype = None
         kwargs = {
            'color'     : eval(self.subsetcolor.get()),
            'nbin'      : eval(self.subsetnbin.get()),
            'fontsize'  : eval(self.subsetfontsize.get()),
            'label'     : eval(self.subsetlabel.get()),
            'ls'        : eval(self.subsetls.get()),
            'marker'    : eval(self.subsetmarker.get()),
            'markersize': eval(self.subsetmarkersize.get()),
            'fillstyle' : eval(self.subsetfillstyle.get()),
            'labelpad'  : eval(self.subsetlabelpad.get()),
         }         
         self.ztfp.add_subset(
            syntax=syntax, astype=astype, **kwargs            
         )
         self.subsetinfo.set('Subset:\n%s' % ( '\n'.join(self.ztfp.dset.keys())))
         
         # update menu
         def updateds(todata, val): todata.set(val)
         # nd
         datasets = list(self.ztfp.dset.keys())
         self.pdataset.set( datasets[0] )
         menu = self.scatterndMenu["menu"]
         menu.delete(0, "end")
         for t in datasets:
            menu.add_command(label=t,command=lambda k=t: updateds(self.pdataset, k))

         # show all
         self.showalldataset.set( datasets[0] )
         menu = self.showallMenu["menu"]
         menu.delete(0, "end")
         for t in datasets:
            menu.add_command(label=t,command=lambda k=t: updateds(self.showalldataset, k))
         
         # table
         self.tabdataset.set( datasets[0] )
         menu = self.stabMenu["menu"]
         menu.delete(0, "end")
         for t in datasets:
            menu.add_command(label=t,command=lambda k=t: updateds(self.tabdataset, k))
            
      inputFrameBut = Tk.Button(self.subsettop, text="Add", command=addsubset)
      inputFrameBut.grid(row=8, column=0, sticky=Tk.NW)
      
      def delsubset():
         if len(self.ztfp.dset) == 0: return        
         self.ztfp.dset.popitem()
         self.ztfp.plotargs.popitem()
         self.subsetinfo.set('Subset:\n%s' % ( '\n'.join(self.ztfp.dset.keys())))

         # update menu
         def updateds(todata, val): self.pdataset.set(todata, val)
         # nd
         datasets = list(self.ztfp.dset.keys())
         if len(datasets) == 0: self.pdataset.set( None )
         else:  self.pdataset.set( datasets[0] )
         menu = self.scatterndMenu["menu"]
         menu.delete(0, "end")
         if len(datasets) > 0:
            for t in datasets:
               menu.add_command(label=t,command=lambda k=t: updateds(self.pdataset, k))

         # show all         
         if len(datasets) == 0: self.showalldataset.set( None )
         else:  self.showalldataset.set( datasets[0] )
         menu = self.showallMenu["menu"]
         menu.delete(0, "end")
         if len(datasets) > 0:
            for t in datasets:
               menu.add_command(label=t,command=lambda k=t: updateds(self.showalldataset, k))
         
         # table         
         if len(datasets) == 0: self.tabdataset.set( None )
         else:  self.tabdataset.set( datasets[0] )
         menu = self.stabMenu["menu"]
         menu.delete(0, "end")
         if len(datasets) > 0:
            for t in datasets:
               menu.add_command(label=t,command=lambda k=t: updateds(self.tabdataset, k))
      inputFrameBut = Tk.Button(self.subsettop, text="Del", command=delsubset)
      inputFrameBut.grid(row=8, column=0, sticky=Tk.N)            
      
      def showmeta():         
         showmeta = self.ztfp.meta
         self.Message( tabulate(showmeta, headers = 'keys', tablefmt = 'grid') )
      inputFrameBut = Tk.Button(self.subsettop, text="Show", command=showmeta)
      inputFrameBut.grid(row=8, column=0, sticky=Tk.NE)
      
   ### popup for Data ###
   def data_popup(self, objid):            
      # pop up window
      self.datatop = Tk.Toplevel(self.parent)
      self.datatop.geometry("500x500+50+50")      
      self.datatop.title('Parse data for %s'%objid)
      
      # data infos
      row1 = Tk.Frame(self.datatop)      
      row1.grid(row=0,column=0,sticky=Tk.W)      
      labelFrame = Tk.LabelFrame(row1, text="Source")
      labelFrame.grid(row=0,column=0)
      
      self.datainfo1 = Tk.StringVar()
      self.datainfo1.set(self.checkdata(objid, ['data']))
      tb = Tk.Label(labelFrame,textvariable=self.datainfo1, fg="red", justify= Tk.LEFT)
      tb.grid(row=0,column=0)
      
      # row1 right
      row11 = Tk.Frame(self.datatop)      
      row11.grid(row=0,column=1,columnspan=2,sticky=Tk.W)
      
      def showauth():
         keypairs = read_default.get_keypairs()
         self.Message( keypairs )
      def showmeta():
         msg = ''
         for k in ['objid', 'z', 'ra', 'dec', 'mkwebv', 'hostebv',
                   'sntype', 'dm', 't0', 'texp', 'fpeak']:
            msg += '%s: %s \n' % (k, str(self.ztfp.data[objid].__dict__[k]))
         self.Message( msg )
      def reset():
         if 'lc' in self.ztfp.data[objid].__dict__:            
            del self.ztfp.data[objid].__dict__['lc']
         if 'spec' in self.ztfp.data[objid].__dict__:            
            del self.ztfp.data[objid].__dict__['spec']
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      authbutton = Tk.Button(row11,text="auth",command=showauth)
      authbutton.grid(row=0, column=0)      
      metabutton = Tk.Button(row11,text="meta",command=showmeta)
      metabutton.grid(row=0, column=1)
      metabutton = Tk.Button(row11,text="reset",command=reset)
      metabutton.grid(row=0, column=2)
      
      if snobject_pars['mjdstart'] is not None:
         mjdstart = '%.1f'%(snobject_pars['mjdstart'])
      elif self.ztfp.data[objid].t0 > 2400000:
         mjdstart = '%.1f' %(float(self.ztfp.data[objid].t0) - 2400000.5 + float(snobject_pars['dstart']))
      else:
         mjdstart = 59000
      if snobject_pars['mjdend'] is not None:
         mjdend = '%.1f'%snobject_pars['mjdend']
      elif self.ztfp.data[objid].t0 > 2400000:
         mjdend = '%.1f'%(float(self.ztfp.data[objid].t0) - 2400000.5 + float(snobject_pars['dend']))
      else:
         mjdend = 59700
      l = Tk.Label(row11,text='MJD:',fg="black")      
      l.grid(row=1,column=0)
      self.minmjdentry = Tk.Entry(row11, justify=Tk.LEFT, width=8)      
      self.minmjdentry.insert(0, mjdstart)
      self.minmjdentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(row11,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.maxmjdentry = Tk.Entry(row11, justify=Tk.LEFT, width=8)
      self.maxmjdentry.insert(0, mjdend)
      self.maxmjdentry.grid(row=1, column=2, sticky=Tk.SW)
      
      snrt = snobject_pars['snrt']
      l = Tk.Label(row11,text='SNR:',fg="black")      
      l.grid(row=2,column=0)
      self.snrentry = Tk.Entry(row11, justify=Tk.LEFT, width=5)
      self.snrentry.insert(0, snrt)
      self.snrentry.grid(row=2, column=1, sticky=Tk.SW)
      l = Tk.Label(row11,text='sigma',fg="black")      
      l.grid(row=2,column=2)
      
      tdbin = snobject_pars['tdbin']
      l = Tk.Label(row11,text='bin:',fg="black")      
      l.grid(row=3,column=0)
      self.tdbentry = Tk.Entry(row11, justify=Tk.LEFT, width=5)
      self.tdbentry.insert(0, tdbin)
      self.tdbentry.grid(row=3, column=1, sticky=Tk.SW)
      l = Tk.Label(row11,text='days',fg="black")      
      l.grid(row=3,column=2)
      
      ### list of options
      # ztf force
      def loadztfforce():
         self.ztfp.data[objid].get_fp_ztf(snrt=float(self.snrentry.get()))
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadztfforce():
         if self.clobberVar.get()==0: clobber = False
         else: clobber = True 
         self.ztfp.data[objid].query_fp_ztf(mjdstart=float(self.minmjdentry.get()),
                                            mjdend=float(self.maxmjdentry.get()), clobber=clobber)
         self.Message('manual input ztf forced photometry file...')         
      data1 = Tk.Label(self.datatop,text='ZTF forced Photometry',fg="black")
      data1.grid(row=1,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadztfforce)
      data1button.grid(row=1, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadztfforce)
      data2button.grid(row=1, column=2)
      # atlas force
      def loadatlasforce():
         if self.clobberVar.get()==0: clobber = False
         else: clobber = True
         if is_number(self.tdbentry.get()): bind = int(self.tdbentry.get())
         else: bind = None
         self.ztfp.data[objid].get_fp_atlas(tdbin=bind, clobber=clobber, snrt=float(self.snrentry.get()))
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadatlasforce():
         if self.clobberVar.get()==0: clobber = False
         else: clobber = True 
         self.ztfp.data[objid].query_fp_atlas(mjdstart=float(self.minmjdentry.get()),
                        mjdend=float(self.maxmjdentry.get()), clobber=clobber)         
      data1 = Tk.Label(self.datatop,text='ATLAS forced Photometry',fg="black")
      data1.grid(row=2,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadatlasforce)
      data1button.grid(row=2, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadatlasforce)
      data2button.grid(row=2, column=2)
      # ztf gm
      def loadztfgm():
         self.ztfp.data[objid].get_alert_ztf(source='marshal', snrt=float(self.snrentry.get()))
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadztfgm():
         self.ztfp.data[objid].query_alert_ztf(source='marshal')                  
      data1 = Tk.Label(self.datatop,text='ZTF Growth marshal Photometry',fg="black")
      data1.grid(row=3,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadztfgm)
      data1button.grid(row=3, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadztfgm)
      data2button.grid(row=3, column=2)
      # ztf fritz
      def loadztffritz():
         self.ztfp.data[objid].get_alert_ztf(source='fritz', snrt=float(self.snrentry.get()))
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadztffritz():
         self.ztfp.data[objid].query_alert_ztf(source='fritz')
      data1 = Tk.Label(self.datatop,text='ZTF fritz Photometry',fg="black")
      data1.grid(row=4,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadztffritz)
      data1button.grid(row=4, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadztffritz)
      data2button.grid(row=4, column=2)
      # open sn catalog photometry
      def loadoac(which):
         self.ztfp.data[objid].get_oac(db='AL', which=which, snrt=float(self.snrentry.get()))
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadoac(which):
         if self.clobberVar.get()==0: clobber = False
         else: clobber = True 
         self.ztfp.data[objid].query_oac(db='AL', which=which, clobber=clobber)
      data1 = Tk.Label(self.datatop,text='Open Astronomy Catalog photometry',fg="black")
      data1.grid(row=5,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=lambda k='photometry': loadoac(k))
      data1button.grid(row=5, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=lambda k='photometry': downloadoac(k))
      data2button.grid(row=5, column=2)
      # GM spectra
      def loadgmspec():
         self.ztfp.data[objid].get_local_spectra(source='marshal')
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadgmspec():
         self.ztfp.data[objid].query_spectra(source='marshal')         
      data1 = Tk.Label(self.datatop,text='ZTF Growth marshal Spectra',fg="black")
      data1.grid(row=6,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadgmspec)
      data1button.grid(row=6, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadgmspec)
      data2button.grid(row=6, column=2)
      # fritz spectra
      def loadfritzspec():
         self.ztfp.data[objid].get_local_spectra(source='fritz')
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")         
      def downloadfritzspec():
         self.ztfp.data[objid].query_spectra(source='fritz')         
      data1 = Tk.Label(self.datatop,text='ZTF fritz Spectra',fg="black")
      data1.grid(row=7,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadfritzspec)
      data1button.grid(row=7, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadfritzspec)
      data2button.grid(row=7, column=2)
      # open sn catalog spectra      
      data1 = Tk.Label(self.datatop,text='Open Astronomy Catalog spectra',fg="black")
      data1.grid(row=8, column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=lambda k='spectra': loadoac(k))
      data1button.grid(row=8, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=lambda k='spectra': downloadoac(k))
      data2button.grid(row=8, column=2)      
      
      # external photometry
      def loadextphot():
         filetypes = [
            ('text files ', '*.txt'),
            ("CSV Files", "*.csv"),
         ]
         filename = filedialog.askopenfilename(parent=self.datatop, initialdir=os.getcwd(),
                                          title='load external photometry', filetypes=filetypes)
         if filename:
            source = self.extphotentry.get()
            if len(source) == 0:
               self.Message('Define a proper source name first')
               return           
            self.ztfp.data[objid].get_external_phot(filename, source, snrt=float(self.snrentry.get()))
            self.datainfo1.set(self.checkdata(objid, ['data']))
            tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")            
      data1 = Tk.Label(self.datatop,text='external Photometry',fg="black")
      data1.grid(row=9,column=0)
      dataFrame = Tk.LabelFrame(self.datatop, text="source name")
      dataFrame.grid(row=9,column=1)
      self.extphotentry = Tk.Entry(dataFrame, width=10)      
      self.extphotentry.grid(row=0, column=0)   
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadextphot)
      data1button.grid(row=9, column=2)         
      # external spectra
      def loadextspec():
         filetypes = [
            ('text files ', '*.txt'),
            ("CSV Files", "*.csv"),
         ]      
         filename = filedialog.askopenfilename(parent=self.datatop, initialdir=os.getcwd(),
                                          title='load external spectra', filetypes=filetypes)            
         if filename:
            epoch = self.extspecentry.get()
            if len(epoch) == 0:
               self.Message('Define a proper epoch first')
               return
            self.ztfp.data[objid].get_external_spectra(filename, epoch, tel='')
            self.datainfo1.set(self.checkdata(objid, ['data']))
            tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")            
      data1 = Tk.Label(self.datatop,text='external Spectra',fg="black")
      data1.grid(row=10,column=0)
      dataFrame = Tk.LabelFrame(self.datatop, text="epochs")
      dataFrame.grid(row=10,column=1)
      self.extspecentry = Tk.Entry(dataFrame, width=10)
      self.extspecentry.grid(row=0, column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadextspec)
      data1button.grid(row=10, column=2)
      # change meta
      def InputMeta():        
         _str = self.inputFrameTxt.get()
         if _str == '': 
            self.Message('Manully input meta for an obj, e.g. z=0.1, sntype=SN Ia, etc')
            return
         if not '=' in _str:
            self.Message('Warning: use "=" to split the key and value, e.g. z=0.1, sntype=SN Ia, etc')
            return
         _key, _value = _str.split('=')
         _key = _key.strip()
         _value = _value.strip()             
         try:
            self.ztfp.data[objid].__dict__[_key] = float(_value)
            self.Message('set %s (type: float) as %s for %s' % (_value, _key, objid))
         except:
            try:
               self.ztfp.data[objid].__dict__[_key] = _value      
               self.Message('set %s (type: str) as %s for %s' % (_value, _key, objid))
            except:
               self.Message('Error: cannot set %s as %s for %s' % (_value, _key, objid))         
      data1 = Tk.Label(self.datatop,text='Change Meta',fg="black")
      data1.grid(row=11,column=0)
      self.inputFrameTxt = Tk.Entry(self.datatop, width=10)
      self.inputFrameTxt.grid(row=11, column=1, sticky=Tk.SW)
      inputFrameBut = Tk.Button(self.datatop,text="submit",command=InputMeta)
      inputFrameBut.grid(row=11, column=2, sticky=Tk.SE)
      
      def midtexp():
         self.ztfp.data[objid].set_texp_midway()         
      inputFrameBut = Tk.Button(self.datatop,text="set texp with (detecitons)",command=midtexp)
      inputFrameBut.grid(row=13, column=0, sticky=Tk.W)

      def cutlc():         
         self.ztfp.data[objid].cut_lc(mjdstart=float(self.minmjdentry.get()),
                                      mjdend=float(self.maxmjdentry.get()))
         self.Message('Cut LC to %.2f - %.2f' %
                      (float(self.minmjdentry.get()), float(self.maxmjdentry.get())))
      inputFrameBut = Tk.Button(self.datatop,text="cut LC",command=cutlc)
      inputFrameBut.grid(row=12, column=2, sticky=Tk.SE)
      
      def checklc():         
         _ne = self.ztfp.data[objid]._nepochs()
         if is_number(self.tdbentry.get()): bind = int(self.tdbentry.get())
         else: bind = 1
         _nc = self.ztfp.data[objid]._ncolors(tdbin=bind)
         _pa = self.ztfp.data[objid]._peak_accuracy(within=3)
         _ep = self.ztfp.data[objid]._earlypoints()
         self.Message('Photometric epochs: %s\nColour epoch: %i\nPeak accuaracy: %s\nEarly points%s' %
                      (_ne, _nc, _pa, _ep))
      inputFrameBut = Tk.Button(self.datatop,text="check LC",command=checklc)
      inputFrameBut.grid(row=12, column=1, sticky=Tk.SW)
      
      dataFrame = Tk.LabelFrame(self.datatop, text="clip sigma")
      dataFrame.grid(row=12,column=0, sticky=Tk.SW)
      self.clipsigma = Tk.Entry(dataFrame, width=10)      
      self.clipsigma.grid(row=0, column=0)      
      def cliplc():
         sigma = None
         if is_number(self.clipsigma.get()): sigma = float(self.clipsigma.get())
         self.ztfp.data[objid].clip_lc(clipsigma=sigma)
      inputFrameBut = Tk.Button(self.datatop,text="clip LC",command=cliplc)
      inputFrameBut.grid(row=12, column=0, sticky=Tk.SE)
      
      def showlc():
         self.Message(tabulate(self.ztfp.data[objid].lc, headers = 'keys', tablefmt = 'grid'), disable=False)
         
      inputFrameBut = Tk.Button(self.datatop,text="LC table",command=showlc)
      inputFrameBut.grid(row=13, column=1, sticky=Tk.SW)
      
      def rapid():         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'rapid'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].rapid(ax=ax)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Rapid plots for %s'%objid)
      inputFrameBut = Tk.Button(self.datatop,text="SN type",command=rapid)
      inputFrameBut.grid(row=13, column=2, sticky=Tk.SE)
      
   ### popup for Data 1 ###
   def data_popup1(self, objid):      
      # pop up window
      self.datatop1 = Tk.Toplevel(self.parent)
      self.datatop1.geometry("500x400+50+50")      
      self.datatop1.title('Fit LC data for %s'%objid)      
      
      # data infos
      row1 = Tk.Frame(self.datatop1)      
      row1.grid(row=0,column=0,sticky=Tk.W)      
      labelFrame = Tk.LabelFrame(row1, text="Fits")
      labelFrame.grid(row=0,column=0)
      
      self.datainfo2 = Tk.StringVar()
      self.datainfo2.set(self.checkdata(objid, ['fits']))
      tb = Tk.Label(labelFrame,textvariable=self.datainfo2, fg="red", justify= Tk.LEFT)
      tb.grid(row=0,column=0)         
      
      # MCMC
      # row1 right
      Frame = Tk.Frame(self.datatop1)
      Frame.grid(row=0,column=1,sticky=Tk.W)
      
      nsteps = snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=1,column=0)
      self.nstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepsentry.insert(0, nsteps)
      self.nstepsentry.grid(row=1, column=1)
      
      nsteps_burnin = snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=2,column=0)
      self.nstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepburnsentry.insert(0, nsteps_burnin)
      self.nstepburnsentry.grid(row=2, column=1)

      nwalkers = snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=1,column=2)
      self.nwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nwalkerentry.insert(0, nwalkers)
      self.nwalkerentry.grid(row=1, column=3)
      
      if 'lc' in self.ztfp.data[objid].__dict__:
         colours = list(np.unique(self.ztfp.data[objid].lc['filter']))
      else:
         colours = snobject_pars['plot_bands']
      l = Tk.Label(Frame,text='colours=',fg="black")      
      l.grid(row=2,column=2)
      self.colorentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.colorentry.insert(0, colours)
      self.colorentry.grid(row=2, column=3)

      # source
      _Frame = Tk.LabelFrame(Frame, text="Source")
      _Frame.grid(row=0,column=0)
      sources = [None]
      if 'lc' in self.ztfp.data[objid].__dict__:
         sources = np.append(sources, np.unique(self.ztfp.data[objid].lc['source']))
      self.lcsource = Tk.StringVar()
      self.lcsource.set( sources[0] )
      
      lcsourceMenu = Tk.OptionMenu(_Frame, self.lcsource, *sources)
      lcsourceMenu.config(justify=Tk.LEFT,width=4)
      lcsourceMenu.grid(row=0,column=0,sticky=Tk.E)
      
      def updatesource(k): self.lcsource.set( k )
      def updateinfos():
         sources = [None]
         if 'lc' in self.ztfp.data[objid].__dict__:            
            sources = np.append(sources, np.unique(self.ztfp.data[objid].lc['source']))
         self.lcsource.set( sources[0] )         
         # update option menu
         menu = lcsourceMenu["menu"]
         menu.delete(0, "end")
         for t in sources:            
            menu.add_command(label=t,command=lambda k=t: updatesource(k))  
      button = Tk.Button(Frame,text="update",command=updateinfos)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=1, sticky=Tk.SW)

      def reset():
         if 'gpcls' in self.ztfp.data[objid].__dict__:
            del self.ztfp.data[objid].__dict__['gpcls']
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'multiband_early' in self.ztfp.data[objid].__dict__['fitcls']:
               del self.ztfp.data[objid].__dict__['fitcls']['multiband_early']
            if 'multiband_main' in self.ztfp.data[objid].__dict__['fitcls']:
               del self.ztfp.data[objid].__dict__['fitcls']['multiband_main']
         self.datainfo2.set(self.checkdata(objid, ['fits']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo2, fg="red", justify= Tk.LEFT)
      button = Tk.Button(Frame,text="reset",command=reset)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=2, sticky=Tk.SW)
      
      # Gaussian process
      Frame = Tk.LabelFrame(self.datatop1, text="Gaussian Process",fg="green")
      Frame.grid(row=1,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['gp_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.gpfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.gpfitminentry.insert(0, gpmin)
      self.gpfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.gpfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.gpfitmaxentry.insert(0, gpmax)
      self.gpfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['gp_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.gpplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.gpplotminentry.insert(0, gpmin)
      self.gpplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.gpplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.gpplotmaxentry.insert(0, gpmax)
      self.gpplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.gproutine = Tk.StringVar()
      self.gproutine.set( snobject_pars['gp_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.gproutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Kernel")
      _Frame.grid(row=1,column=3)
      _list = ['matern52', 'matern32', 'squaredexp']
      self.kernel = Tk.StringVar()
      self.kernel.set( snobject_pars['kernel'] )
      Menu = Tk.OptionMenu(_Frame, self.kernel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Mean")
      _Frame.grid(row=0,column=4)
      _list = ['mean', 'gaussian', 'bazin', 'villar']
      self.gpmean = Tk.StringVar()
      self.gpmean.set( snobject_pars['gp_mean'] )
      Menu = Tk.OptionMenu(_Frame, self.gpmean, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)   
      
      _Frame = Tk.LabelFrame(Frame, text="Fix scale")
      _Frame.grid(row=1,column=4)
      _list = [False, True]
      self.gpfixscale = Tk.StringVar()
      self.gpfixscale.set( snobject_pars['fix_scale'] )
      Menu = Tk.OptionMenu(_Frame, self.gpfixscale, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)   
      
      # Power law fits
      Frame = Tk.LabelFrame(self.datatop1, text="Multiband early",fg="blue")
      Frame.grid(row=2,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['multiband_early_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.plfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.plfitminentry.insert(0, gpmin)
      self.plfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.plfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.plfitmaxentry.insert(0, gpmax)
      self.plfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['multiband_early_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.plplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.plplotminentry.insert(0, gpmin)
      self.plplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.plplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.plplotmaxentry.insert(0, gpmax)
      self.plplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.plroutine = Tk.StringVar()
      self.plroutine.set( snobject_pars['multiband_early_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.plroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      #_list = ['pl', 'pl2']
      _list = get_model(which_engine='multiband_early', with_alias=False)['multiband_early']
      self.plmodel = Tk.StringVar()
      self.plmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.plmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # Bazin fits
      Frame = Tk.LabelFrame(self.datatop1, text="Multiband main",fg="brown")
      Frame.grid(row=3,column=0,columnspan=2,sticky=Tk.W)

      gpmin, gpmax = snobject_pars['multiband_main_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.bzfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bzfitminentry.insert(0, gpmin)
      self.bzfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.bzfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bzfitmaxentry.insert(0, gpmax)
      self.bzfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['multiband_main_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.bzplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bzplotminentry.insert(0, gpmin)
      self.bzplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.bzplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bzplotmaxentry.insert(0, gpmax)
      self.bzplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.bzroutine = Tk.StringVar()
      self.bzroutine.set( snobject_pars['multiband_main_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.bzroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      #_list = ['bz', 'bz2']
      _list = get_model(which_engine='multiband_main', with_alias=False)['multiband_main']
      self.bzmodel = Tk.StringVar()
      self.bzmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.bzmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # done
      def rungp():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.lcsource.get()
         self.Message('GP for source: %s with clobber=%s'%(source,clobber))
         if len(source) == 0 or source == 'None': source=None
         gp_bands = self.colorentry.get().split()
         gp_fitr = [float(self.gpfitminentry.get()),float(self.gpfitmaxentry.get())]
         gp_plotr = [float(self.gpplotminentry.get()),float(self.gpplotmaxentry.get())]
         self.ztfp.data[objid].run_gp(
            gp_fit = True,
            gp_bands = gp_bands,
            gp_xrange = gp_fitr,
            gp_xrangep = gp_plotr,
            kernel = self.kernel.get(),
            fix_scale = self.gpfixscale.get(),
            gp_mean = self.gpmean.get(),
            gp_routine = self.gproutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            gp_redo = clobber,
            source = source,
         )
         self.datainfo2.set(self.checkdata(objid, ['fits']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo2, fg="red", justify= Tk.LEFT)
      button = Tk.Button(self.datatop1,text="run",command=rungp,fg="green")
      button.grid(row=1, column=2, sticky=Tk.NE)
      def runpl():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.lcsource.get()
         self.Message('Multi early fits for source: %s with clobber=%s'%(source,clobber))
         if len(source) == 0 or source == 'None': source=None
         fit_bands = self.colorentry.get().split()
         fit_fitr = [float(self.plfitminentry.get()),float(self.plfitmaxentry.get())]
         fit_fitp = [float(self.plplotminentry.get()),float(self.plplotmaxentry.get())]   
         self.ztfp.data[objid].run_fit(
            'multiband_early',
            fit_methods = [self.plmodel.get()],
            multiband_early_bands = fit_bands,
            multiband_early_xrange = fit_fitr,
            multiband_early_xrangep = fit_fitp,
            multiband_early_routine = self.plroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            fit_redo = clobber,
            source = source,
         )
         self.datainfo2.set(self.checkdata(objid, ['fits']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo2, fg="red", justify= Tk.LEFT)
      button1 = Tk.Button(self.datatop1,text="run",command=runpl,fg="blue")
      button1.grid(row=2, column=1, sticky=Tk.NE)
      def runbz():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.lcsource.get()
         self.Message('Multi main fits for source: %s with clobber=%s'%(source,clobber))
         if len(source) == 0 or source == 'None': source=None
         fit_bands = self.colorentry.get().split()
         fit_fitr = [float(self.bzfitminentry.get()),float(self.bzfitmaxentry.get())]
         fit_fitp = [float(self.bzplotminentry.get()),float(self.bzplotmaxentry.get())]
         self.ztfp.data[objid].run_fit(
            'multiband_main',
            fit_methods = [self.bzmodel.get()],
            multiband_main_bands = fit_bands,
            multiband_main_xrange = fit_fitr,
            multiband_main_xrangep = fit_fitp,
            multiband_main_routine = self.bzroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            fit_redo = clobber,
            source = source,
         )
         self.datainfo2.set(self.checkdata(objid, ['fits']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo2, fg="red", justify= Tk.LEFT)         
      button2 = Tk.Button(self.datatop1,text="run",command=runbz,fg="brown")
      button2.grid(row=3, column=1, sticky=Tk.NE)
      
      #
      Frame = Tk.LabelFrame(self.datatop1, text="set t")
      Frame.grid(row=2,column=2,rowspan=2, sticky=Tk.W)
      
      def gpt0(): self.ztfp.data[objid].set_peak_gp(snobject_pars['set_tpeak_filter'])
      inputFrameBut = Tk.Button(Frame, text="GP: t0", command=gpt0,fg="green")
      inputFrameBut.grid(row=0, column=0, sticky=Tk.NW)

      def maint0():
         self.ztfp.data[objid].set_peak_multiband_main(snobject_pars['set_tpeak_filter'])
      inputFrameBut = Tk.Button(Frame, text="main: t0", command=maint0,fg="blue")
      inputFrameBut.grid(row=1, column=0, sticky=Tk.W)      

      def pltexp():
         self.ztfp.data[objid].set_texp_pl(snobject_pars['set_texp_filter'])
      inputFrameBut = Tk.Button(Frame, text="early: texp", command=pltexp,fg="brown")
      inputFrameBut.grid(row=2, column=0, sticky=Tk.SW)
      
      # corner plots
      def gpcontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, gp=True, clobber=clobber)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)         
      inputFrameBut = Tk.Button(self.datatop1, text="Contours", command=gpcontours,fg="green")
      inputFrameBut.grid(row=1, column=2, sticky=Tk.SE)

      def maincontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True         
         fit_bands = self.colorentry.get().split()
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass         
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')         
         self.ztfp.data[objid].show_corner(ax, engine='multiband_main',
                           model=self.bzmodel.get(), clobber=clobber, filts=fit_bands)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.datatop1, text="Contours", command=maincontours,fg="brown")
      inputFrameBut.grid(row=3, column=1, sticky=Tk.SE)
      
      def plcontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         fit_bands = self.colorentry.get().split()
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, engine='multiband_early',
                           model=self.plmodel.get(), clobber=clobber, filts=fit_bands)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim()) 
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.datatop1, text="Contours", command=plcontours,fg="blue")
      inputFrameBut.grid(row=2, column=1, sticky=Tk.SE)      

      # multiband early description
      def medescription():
         model=self.plmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.datatop1, text="describe", command=medescription, fg='blue')
      inputFrameBut.grid(row=2, column=1, sticky=Tk.E)
      
      # multiband main description
      def mmdescription():
         model=self.bzmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.datatop1, text="describe", command=mmdescription, fg='brown')
      inputFrameBut.grid(row=3, column=1, sticky=Tk.E)
      
   ### popup for bb  ###
   def bb_popup(self, objid):
            
      # pop up window
      self.bbtop = Tk.Toplevel(self.parent)
      self.bbtop.geometry("600x600+50+50")      
      self.bbtop.title('Make color/bolometric for %s'%objid)
      
      # data infos
      row1 = Tk.Frame(self.bbtop)      
      row1.grid(row=0,column=0,sticky=Tk.W)      
      labelFrame = Tk.LabelFrame(row1, text="Data")
      labelFrame.grid(row=0,column=0)
      
      self.datainfo3 = Tk.StringVar()
      self.datainfo3.set(self.checkdata(objid, ['bc']))
      tb = Tk.Label(labelFrame,textvariable=self.datainfo3,fg="red")      
      tb.grid(row=0,column=0)

      # row1 middle
      row11 = Tk.Frame(self.bbtop)
      row11.grid(row=0,column=1,sticky=Tk.W)      
                  
      row = Tk.LabelFrame(row11, text="Match\ncolours")
      row.grid(row=0,column=0)
      
      self.interplist = Tk.Listbox(row, selectmode = "multiple")  
      self.interplist.grid(row=0,column=0)
      self.interplist.config(width=8, height=3)
      for each_item in ["bin","gp","fit"]: self.interplist.insert(Tk.END, each_item)
      self.interplist.selection_set( first = 0 )
      
      # row1 right
      row11 = Tk.Frame(self.bbtop)
      row11.grid(row=0,column=2,sticky=Tk.W)

      snrt = snobject_pars['snrt']
      l = Tk.Label(row11,text='SNR:',fg="black")      
      l.grid(row=0,column=0)
      self.snrentry1 = Tk.Entry(row11, justify=Tk.LEFT, width=5)
      self.snrentry1.insert(0, snrt)
      self.snrentry1.grid(row=0, column=1, sticky=Tk.SW)
      l = Tk.Label(row11,text='sigma',fg="black")      
      l.grid(row=0,column=2)
      
      tdbin = snobject_pars['tdbin']
      l = Tk.Label(row11,text='bin:',fg="black")      
      l.grid(row=1,column=0)
      self.tdbentry1 = Tk.Entry(row11, justify=Tk.LEFT, width=5)
      self.tdbentry1.insert(0, tdbin)
      self.tdbentry1.grid(row=1, column=1, sticky=Tk.SW)
      l = Tk.Label(row11,text='days',fg="black")      
      l.grid(row=1,column=2)

      self.corrmkwebv = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(row11, text="corr mkw",variable=self.corrmkwebv) 
      cButton.grid(row=2,column=0,columnspan=3,sticky=Tk.W)
      
      self.corrhostebv = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(row11, text="corr host",variable=self.corrhostebv)      
      cButton.grid(row=3,column=0,columnspan=3,sticky=Tk.W)

      def reset():         
         if 'colors' in self.ztfp.data[objid].__dict__:
            del self.ztfp.data[objid].__dict__['colors']
         if 'mbol' in self.ztfp.data[objid].__dict__:            
            del self.ztfp.data[objid].__dict__['mbol']
         if 'mbolbb' in self.ztfp.data[objid].__dict__:            
            del self.ztfp.data[objid].__dict__['mbolbb']
         if 'mbolspec' in self.ztfp.data[objid].__dict__:            
            del self.ztfp.data[objid].__dict__['mbolspec']
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo3,fg="red")                
      button = Tk.Button(row11,text="reset",command=reset)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=4,column=0,columnspan=3,sticky=Tk.W)
      
      # calculate colours
      if snobject_pars['color_bands'] is None:
         colours = ['g','r']
      else:
         colours = snobject_pars['color_bands']              
      data1 = Tk.Label(self.bbtop,text='calculate colors',fg="black")
      data1.grid(row=1,column=0)        
      colourFrame = Tk.LabelFrame(self.bbtop, text="2 bands")
      colourFrame.grid(row=1,column=1)
      self.colorentry1 = Tk.Entry(colourFrame, width=10)
      self.colorentry1.insert(0, colours)
      self.colorentry1.grid(row=0, column=0)      
      def makecolor():
         bands = self.colorentry1.get().split()                           
         if len(bands) != 2:
            self.Message('Need/Only two bands needed')
            return
         cinterp = []
         for index in self.interplist.curselection():
            cinterp.append(self.interplist.get(index))   
         if len(cinterp) == 0:
            self.Message('Select at least one option to match colours')
            return
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         self.ztfp.data[objid].calc_colors(
            color_bands = bands,
            tdbin = int(self.tdbentry1.get()),
            color_interp = cinterp,
            snrt = int(self.snrentry1.get()),
            corr_mkw = corr_mkw,
            corr_host = corr_host,
         )
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(self.bbtop,textvariable=self.datainfo3,fg="red")
      data1button = Tk.Button(self.bbtop,text="submit",command=makecolor)
      data1button.grid(row=1, column=2)
      
      # make lyman bolometric
      if snobject_pars['color_bands'] is None:
         colours = ['g','r']
      else:
         colours = snobject_pars['color_bands'] 
      data1 = Tk.Label(self.bbtop,text='Bolometric with Lyman BC',fg="black")
      data1.grid(row=2,column=0)
      lymanFrame = Tk.LabelFrame(self.bbtop, text="2 bands")
      lymanFrame.grid(row=2,column=1)      
      self.lymanentry = Tk.Entry(lymanFrame, width=10)
      self.lymanentry.insert(0, colours)
      self.lymanentry.grid(row=0, column=0)
      
      def makebc():                           
         bands = self.lymanentry.get().split()
         if len(bands) != 2:
            self.Message('Need/Only two bands needed')            
            return
         cinterp = []
         for index in self.interplist.curselection():
            cinterp.append(self.interplist.get(index))
         if len(cinterp) == 0:
            self.Message('Select at least one option to match colours')
            return
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         self.ztfp.data[objid].lyman_bol(
            make_bol = ['lyman'],
            lyman_bands = bands,
            lyman_interp = cinterp,
            tdbin = int(self.tdbentry1.get()),            
            snrt = int(self.snrentry1.get()),
            corr_mkw = corr_mkw,
            corr_host = corr_host,            
         )
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(self.bbtop,textvariable=self.datainfo3,fg="red")         
      data1button = Tk.Button(self.bbtop,text="submit",command=makebc)
      data1button.grid(row=2, column=2)
      
      # make BB bolometric
      if snobject_pars['sed_bands'] is None:
         colours = list(np.unique(self.ztfp.data[objid].lc['filter']))
      else:
         colours = snobject_pars['sed_bands']
      data1 = Tk.Label(self.bbtop,text='Bolometric with BB',fg="black")
      data1.grid(row=3,column=0)        
      bbFrame = Tk.LabelFrame(self.bbtop, text="multi bands")
      bbFrame.grid(row=3,column=1)      
      self.bbentry = Tk.Entry(bbFrame, width=10)
      self.bbentry.insert(0, colours)
      self.bbentry.grid(row=0, column=0)
      
      def makebc1():
         bands = self.bbentry.get()
         if len(bands) < 2:
            self.Message('Need at least three bands')
            return
         cinterp = []
         for index in self.interplist.curselection():
            cinterp.append(self.interplist.get(index))
         if len(cinterp) == 0:
            self.Message('Select at least one option to match colours')
            return
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         self.ztfp.data[objid].bb_colors(
            sed_bands=bands,
            sed_color_interp=cinterp,
            tdbin = int(self.tdbentry1.get()),            
            snrt = int(self.snrentry1.get()),
            corr_mkw = corr_mkw,
            corr_host = corr_host, 
         ) 
         self.ztfp.data[objid].bb_bol(
            fastsedfitting=True,
            make_bol=['bb','spec'],
            sed_bands=bands,
            sed_color_interp=cinterp,
         )
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(self.bbtop,textvariable=self.datainfo3,fg="red") 
      data1button = Tk.Button(self.bbtop,text="submit",command=makebc1)
      data1button.grid(row=3, column=2)
      
      def sedpop():
         if not 'sedtop' in self.__dict__: self.sed_popup(objid)
         elif self.sedtop.winfo_exists() == 0: self.sed_popup(objid)
         else: self.Message('SED fitting topup window already exists')
      data1button = Tk.Button(self.bbtop,text="SED fittings",command=sedpop)
      data1button.grid(row=3, column=3)
      
      # estimate host
      if snobject_pars['hostebv_bands'] is None:
         colours = ['g','r']
      else:
         colours = snobject_pars['hostebv_bands']   
      data1 = Tk.Label(self.bbtop,text='Estimate host ebv',fg="black")
      data1.grid(row=4,column=0)
      hostFrame = Tk.LabelFrame(self.bbtop, text="2 bands")
      hostFrame.grid(row=4,column=1)      
      self.hostentry = Tk.Entry(hostFrame, width=10)
      self.hostentry.insert(0, colours)
      self.hostentry.grid(row=0, column=0)
      
      def esthost():
         bands = self.hostentry.get().split()         
         cinterp = self.interplist.get(Tk.ACTIVE) 
         if len(bands) != 2:
            self.Message('Need/Only two bands needed')
            return
         hostebv = self.ztfp.data[objid].hostebv
         self.ztfp.data[objid].est_hostebv_with_colours(
            index=0, returnv=False,
            hostebv_interp=cinterp,            
            hostebv_bands=bands,
            tdbin = int(self.tdbentry1.get()),            
            snrt = int(self.snrentry1.get()),
         )         
         self.Message( 'With %s, Host EBV: %.2f -> %.2f' % (cinterp, hostebv, self.ztfp.data[objid].hostebv) )
      data1button = Tk.Button(self.bbtop,text="submit",command=esthost)
      data1button.grid(row=4, column=2)

      # mag at
      data1 = Tk.Label(self.bbtop,text='app Mag',fg="black")
      data1.grid(row=5,column=0)
      
      bFrame = Tk.LabelFrame(self.bbtop, text="band")
      bFrame.grid(row=5,column=1)      
      self.appmagentry = Tk.Entry(bFrame, width=10)
      self.appmagentry.insert(0, 'r')
      self.appmagentry.grid(row=0, column=0, sticky=Tk.W)

      bFrame = Tk.LabelFrame(self.bbtop, text="phase")
      bFrame.grid(row=5,column=2)
      self.appmagpentry = Tk.Entry(bFrame, width=10)
      self.appmagpentry.insert(0, 0)
      self.appmagpentry.grid(row=0, column=0, sticky=Tk.E)
      
      def magat():
         filt = self.appmagentry.get()
         phase = float(self.appmagpentry.get())
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         cinterp = self.interplist.get(Tk.ACTIVE)         
         m, em = self.ztfp.data[objid]._mag_at(
            filt, phase,
            interpolation=cinterp,
            corr_mkw=corr_mkw,
            corr_host=corr_host,
         )
         self.Message( 'With %s, Mag = %s (%s)' % (cinterp, m, em) )
      data1button = Tk.Button(self.bbtop,text="submit",command=magat)
      data1button.grid(row=5, column=3)

      # abs mag at        
      data1 = Tk.Label(self.bbtop,text='abs Mag',fg="black")
      data1.grid(row=6,column=0)
      
      bFrame = Tk.LabelFrame(self.bbtop, text="band")
      bFrame.grid(row=6,column=1)      
      self.absmagentry = Tk.Entry(bFrame, width=10)
      self.absmagentry.insert(0, 'r')
      self.absmagentry.grid(row=0, column=0, sticky=Tk.W)

      bFrame = Tk.LabelFrame(self.bbtop, text="phase")
      bFrame.grid(row=6,column=2)
      self.absmagpentry = Tk.Entry(bFrame, width=10)
      self.absmagpentry.insert(0, 0)
      self.absmagpentry.grid(row=0, column=0, sticky=Tk.E)
      
      def absmagat():
         filt = self.absmagentry.get()
         phase = float(self.absmagpentry.get())
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         cinterp = self.interplist.get(Tk.ACTIVE)
         m, em = self.ztfp.data[objid]._absmag_at(
            filt, phase,
            interpolation=cinterp,
            corr_mkw=corr_mkw,
            corr_host=corr_host,
         )
         self.Message( 'With %s, Mag = %s (%s)' % (cinterp, m, em) )
      data1button = Tk.Button(self.bbtop,text="submit",command=absmagat)
      data1button.grid(row=6, column=3)
      
      # delta mag at        
      data1 = Tk.Label(self.bbtop,text='delta Mag p1-p2',fg="black")
      data1.grid(row=7,column=0)
      
      bFrame = Tk.LabelFrame(self.bbtop, text="band")
      bFrame.grid(row=7,column=1)      
      self.dmatentry = Tk.Entry(bFrame, width=10)
      self.dmatentry.insert(0, 'r')
      self.dmatentry.grid(row=0, column=0, sticky=Tk.W)

      bFrame = Tk.LabelFrame(self.bbtop, text="phase")
      bFrame.grid(row=7,column=2)
      self.dmatpentry = Tk.Entry(bFrame, width=10)
      self.dmatpentry.insert(0, [0,15])
      self.dmatpentry.grid(row=0, column=0, sticky=Tk.E)
      
      def dmagat():
         filt = self.dmatentry.get()
         phases = self.dmatpentry.get().split()
         if len(phases) != 2:
            self.Message('Need/Only two phases needed')
            return
         cinterp = self.interplist.get(Tk.ACTIVE)
         m, em = self.ztfp.data[objid]._rate_at(
            filt, float(phases[0]), float(phases[1]),
            interpolation=cinterp,
         )
         self.Message( 'With %s, delta Mag (%s-%s) = %s (%s)' % (cinterp, phases[0], phases[1], m, em) )
      data1button = Tk.Button(self.bbtop,text="submit",command=dmagat)
      data1button.grid(row=7, column=3)

      # delta mag at
      data1 = Tk.Label(self.bbtop,text='delta Mag c1-c2',fg="black")
      data1.grid(row=8,column=0)
      
      bFrame = Tk.LabelFrame(self.bbtop, text="2 bands")
      bFrame.grid(row=8,column=1)      
      self.dmatentry1 = Tk.Entry(bFrame, width=10)
      self.dmatentry1.insert(0, ['g', 'r'])
      self.dmatentry1.grid(row=0, column=0, sticky=Tk.W)
      
      bFrame = Tk.LabelFrame(self.bbtop, text="phase")
      bFrame.grid(row=8,column=2)
      self.dmatpentry1 = Tk.Entry(bFrame, width=10)
      self.dmatpentry1.insert(0, 0)
      self.dmatpentry1.grid(row=0, column=0, sticky=Tk.E)
      
      def cmagat():
         filters = self.dmatentry1.get().split()
         phase = self.dmatpentry1.get()
         if len(filters) != 2:
            self.Message('Need/Only two filters')
            return
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         cinterp=self.interplist.get(Tk.ACTIVE)        
         m, em = self.ztfp.data[objid]._color_at(
            filters[0], filters[1], float(phase),
            interpolation=cinterp,
            corr_mkw=corr_mkw,
            corr_host=corr_host,            
         )
         self.Message( 'With %s, delta Mag (%s-%s) = %s (%s)' % (cinterp, filters[0], filters[1], m, em) )      
      data1button = Tk.Button(self.bbtop,text="submit",command=cmagat)
      data1button.grid(row=8, column=3)
      
   ### popup for bb 1 ###
   def bb_popup1(self, objid):
      
      # pop up window
      self.bbtop1 = Tk.Toplevel(self.parent)
      self.bbtop1.geometry("450x600+50+50")      
      self.bbtop1.title('Bolometric LC fit for %s'%objid)
      
      # data infos
      row1 = Tk.Frame(self.bbtop1)      
      row1.grid(row=0,column=0,sticky=Tk.W)      
      labelFrame = Tk.LabelFrame(row1, text="Bolometric Fits")
      labelFrame.grid(row=0,column=0)
      
      self.datainfo4 = Tk.StringVar()
      self.datainfo4.set(self.checkdata(objid, ['bolfit']))
      tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")      
      tb.grid(row=0,column=0)
      
      # row1 right
      Frame = Tk.Frame(self.bbtop1)
      Frame.grid(row=0,column=1,sticky=Tk.W)
      
      nsteps = snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=1,column=0)
      self.nstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepsentry.insert(0, nsteps)
      self.nstepsentry.grid(row=1, column=1)
      
      nsteps_burnin = snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=2,column=0)
      self.nstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepburnsentry.insert(0, nsteps_burnin)
      self.nstepburnsentry.grid(row=2, column=1)

      nwalkers = snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=1,column=2)
      self.nwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nwalkerentry.insert(0, nwalkers)
      self.nwalkerentry.grid(row=1, column=3)
      
      if snobject_pars['color_interp'] is None:         
         cinterp = ['bin','gp','fit']
      else:
         cinterp = []         
         for _ in snobject_pars['color_interp']:
            cinterp.append(_)
      l = Tk.Label(Frame,text='inter=',fg="black")      
      l.grid(row=2,column=2)
      self.cinterpentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.cinterpentry.insert(0, cinterp)
      self.cinterpentry.grid(row=2, column=3)

      # source
      _Frame = Tk.LabelFrame(Frame, text="Source")
      _Frame.grid(row=0,column=0)
      sources = [None]
      if 'mbol' in self.ztfp.data[objid].__dict__: sources = np.append(sources, 'mbol')
      if 'mbolbb' in self.ztfp.data[objid].__dict__: sources = np.append(sources, 'mbolbb')
      if 'mbolspec' in self.ztfp.data[objid].__dict__: sources = np.append(sources, 'mbolspec')
      self.msource = Tk.StringVar()
      self.msource.set( sources[0] )
      
      lcsourceMenu = Tk.OptionMenu(_Frame, self.msource, *sources)
      lcsourceMenu.config(justify=Tk.LEFT,width=4)
      lcsourceMenu.grid(row=0,column=0,sticky=Tk.E)
      
      def updatesource(k): self.msource.set( k )
      def updateinfos():
         sources = [None]
         if 'mbol' in self.ztfp.data[objid].__dict__: sources = np.append(sources, 'mbol')
         if 'mbolbb' in self.ztfp.data[objid].__dict__: sources = np.append(sources, 'mbolbb')
         if 'mbolspec' in self.ztfp.data[objid].__dict__: sources = np.append(sources, 'mbolspec')
         self.msource.set( sources[0] )         
         # update option menu
         menu = lcsourceMenu["menu"]
         menu.delete(0, "end")
         for t in sources:            
            menu.add_command(label=t,command=lambda k=t: updatesource(k))  
      button = Tk.Button(Frame,text="update",command=updateinfos)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=1, sticky=Tk.SW)

      def reset():         
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            for k in ['bol_early', 'bol_main', 'bol_tail', 'bol_full']:
               if k in self.ztfp.data[objid].__dict__['fitcls']:
                  del self.ztfp.data[objid].__dict__['fitcls'][k]           
         self.datainfo4.set(self.checkdata(objid, ['bolfit']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")   
      button = Tk.Button(Frame,text="reset",command=reset)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=2, sticky=Tk.SW)
      
      # bol early fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric early (shock cooling)")
      Frame.grid(row=1,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_early_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.sbofitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sbofitminentry.insert(0, gpmin)
      self.sbofitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.sbofitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sbofitmaxentry.insert(0, gpmax)
      self.sbofitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_early_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.sboplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sboplotminentry.insert(0, gpmin)
      self.sboplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.sboplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sboplotmaxentry.insert(0, gpmax)
      self.sboplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.sboroutine = Tk.StringVar()
      self.sboroutine.set( snobject_pars['bol_early_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.sboroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      #_list = ['shock']
      _list = get_model(which_engine='bol_early', with_alias=False)['bol_early']
      self.sbomodel = Tk.StringVar()
      self.sbomodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.sbomodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      # bol main fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric main (Arnett)",fg='green')
      Frame.grid(row=2,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_main_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.bmfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bmfitminentry.insert(0, gpmin)
      self.bmfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.bmfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bmfitmaxentry.insert(0, gpmax)
      self.bmfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_main_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.bmplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bmplotminentry.insert(0, gpmin)
      self.bmplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.bmplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bmplotmaxentry.insert(0, gpmax)
      self.bmplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.bmroutine = Tk.StringVar()
      self.bmroutine.set( snobject_pars['bol_main_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.bmroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = get_model(which_engine='bol_main', with_alias=False)['bol_main']
      self.bmmodel = Tk.StringVar()
      #self.bmmodel.set( _list[0] )
      self.bmmodel.set( 'arnett_fit_taum_texp' )
      Menu = Tk.OptionMenu(_Frame, self.bmmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      # bol late fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric tail (radioactive)",fg='blue')
      Frame.grid(row=3,column=0,columnspan=2,sticky=Tk.W)

      gpmin, gpmax = snobject_pars['bol_tail_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.tailfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.tailfitminentry.insert(0, gpmin)
      self.tailfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.tailfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.tailfitmaxentry.insert(0, gpmax)
      self.tailfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_tail_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.tailplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.tailplotminentry.insert(0, gpmin)
      self.tailplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.tailplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.tailplotmaxentry.insert(0, gpmax)
      self.tailplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.tailroutine = Tk.StringVar()
      self.tailroutine.set( snobject_pars['bol_tail_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.tailroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      #_list = ['tail', 'tail2', 'tail3', 'tail4']
      _list = get_model(which_engine='bol_tail', with_alias=False)['bol_tail']
      self.tailmodel = Tk.StringVar()
      self.tailmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.tailmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # bol joint fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric full (Arnett+tail)",fg='brown')
      Frame.grid(row=4,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_full_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.bffitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bffitminentry.insert(0, gpmin)
      self.bffitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.bffitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bffitmaxentry.insert(0, gpmax)
      self.bffitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['bol_full_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.bfplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bfplotminentry.insert(0, gpmin)
      self.bfplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.bfplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.bfplotmaxentry.insert(0, gpmax)
      self.bfplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.bfroutine = Tk.StringVar()
      self.bfroutine.set( snobject_pars['bol_full_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.bfroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)      
      _list = get_model(which_engine='bol_full', with_alias=False)['bol_full']
      self.bfmodel = Tk.StringVar()
      self.bfmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.bfmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # done
      def runearlybol():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.msource.get()
         self.Message('Bol early fits for source: %s'%(source))
         if len(source) == 0 or source == 'None': source=None  
         fit_fitr = [float(self.sbofitminentry.get()),float(self.sbofitmaxentry.get())]
         fit_fitp = [float(self.sboplotminentry.get()),float(self.sboplotmaxentry.get())]
         cinterp = []
         for c in self.cinterpentry.get().split(): cinterp.append(c.strip())
         if len(cinterp) == 0:
            self.Message('Error: select reasoble interpolations: gp/bin/fit')
            return   
         self.ztfp.data[objid].run_fit(
            'bol_early',
            fit_methods = [self.sbomodel.get()],
            bol_early_xrange = fit_fitr,
            bol_early_xrangep = fit_fitp,
            bol_early_routine = self.sboroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            bol_early_color_interp = cinterp,
            fit_redo = clobber,
            source = source,
         )
         self.datainfo4.set(self.checkdata(objid, ['bolfit']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")   
      button = Tk.Button(self.bbtop1,text="Run",command=runearlybol)
      button.grid(row=1, column=1, sticky=Tk.NE)
      
      def runmainbol():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.msource.get()
         self.Message('Bol main fits for source: %s'%(source))
         if len(source) == 0 or source == 'None': source=None  
         fit_fitr = [float(self.bmfitminentry.get()),float(self.bmfitmaxentry.get())]
         fit_fitp = [float(self.bmplotminentry.get()),float(self.bmplotmaxentry.get())]
         cinterp = []
         for c in self.cinterpentry.get().split(): cinterp.append(c.strip())
         if len(cinterp) == 0:
            self.Message('Error: select reasoble interpolations: gp/bin/fit')
            return
         self.ztfp.data[objid].run_fit(
            'bol_main',
            fit_methods = [self.bmmodel.get()],
            bol_main_xrange = fit_fitr,
            bol_main_xrangep = fit_fitp,
            bol_main_routine = self.bmroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            bol_main_color_interp = cinterp,
            fit_redo = clobber,
            source = source,
         )
         self.datainfo4.set(self.checkdata(objid, ['bolfit']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")      
      button1 = Tk.Button(self.bbtop1,text="Run",command=runmainbol,fg='green')
      button1.grid(row=2, column=1, sticky=Tk.NE)

      def runtailbol():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.msource.get()
         self.Message('Bol tail fits for source: %s'%(source))
         if len(source) == 0 or source == 'None': source=None  
         fit_fitr = [float(self.tailfitminentry.get()),float(self.tailfitmaxentry.get())]
         fit_fitp = [float(self.tailplotminentry.get()),float(self.tailplotmaxentry.get())]
         cinterp = []
         for c in self.cinterpentry.get().split(): cinterp.append(c.strip())
         if len(cinterp) == 0:
            self.Message('Error: select reasoble interpolations: gp/bin/fit')
            return          
         self.ztfp.data[objid].run_fit(
            'bol_tail',
            fit_methods = [self.tailmodel.get()],
            bol_tail_xrange = fit_fitr,
            bol_tail_xrangep = fit_fitp,
            bol_tail_routine = self.tailroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            bol_tail_color_interp = cinterp,
            fit_redo=clobber,
            source=source,
         )
         self.datainfo4.set(self.checkdata(objid, ['bolfit']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")      
      button1 = Tk.Button(self.bbtop1,text="Run",command=runtailbol,fg='blue')
      button1.grid(row=3, column=1, sticky=Tk.NE)

      def runfullbol():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.msource.get()
         self.Message('Bol joint fits for source: %s'%(source))
         if len(source) == 0 or source == 'None': source=None  
         fit_fitr = [float(self.bffitminentry.get()),float(self.bffitmaxentry.get())]
         fit_fitp = [float(self.bfplotminentry.get()),float(self.bfplotmaxentry.get())]
         cinterp = []
         for c in self.cinterpentry.get().split(): cinterp.append(c)            
         if len(cinterp) == 0:
            self.Message('Error: select reasoble interpolations: gp/bin/fit')
            return          
         self.ztfp.data[objid].run_fit(
            'bol_full',
            fit_methods = [self.bfmodel.get()],
            bol_tail_xrange = fit_fitr,
            bol_tail_xrangep = fit_fitp,
            bol_tail_routine = self.bfroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            bol_tail_color_interp = cinterp,
            fit_redo=clobber,
            source=source,
         )
         self.datainfo4.set(self.checkdata(objid, ['bolfit']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")      
      button1 = Tk.Button(self.bbtop1,text="Run",command=runfullbol,fg='brown')
      button1.grid(row=4, column=1, sticky=Tk.NE)

      # bol early description
      def bearlydescription():
         model=self.sbomodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.bbtop1, text="describe", command=bearlydescription)
      inputFrameBut.grid(row=1, column=1, sticky=Tk.E)
      
      # bol early corner
      def bearlycontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         source=self.msource.get()
         if len(source) == 0 or source == 'None': source=None           
         self.ztfp.data[objid].show_corner(
            ax, engine='bol_early', model=self.sbomodel.get(), source=source,
            clobber=clobber,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.bbtop1, text="Contour", command=bearlycontours)
      inputFrameBut.grid(row=1, column=1, sticky=Tk.SE)

      # bol main description
      def bmaindescription():
         model=self.bmmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.bbtop1, text="describe", command=bmaindescription, fg='green')
      inputFrameBut.grid(row=2, column=1, sticky=Tk.E)
      
      # bol main corner
      def bmaincontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         source=self.msource.get()
         if len(source) == 0 or source == 'None': source=None           
         self.ztfp.data[objid].show_corner(
            ax, engine='bol_main', model=self.bmmodel.get(), source=source,
            clobber=clobber,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.bbtop1, text="Contour", command=bmaincontours,fg='green')
      inputFrameBut.grid(row=2, column=1, sticky=Tk.SE)

      # bol main description
      def btaildescription():
         model=self.tailmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.bbtop1, text="describe", command=btaildescription, fg='blue')
      inputFrameBut.grid(row=3, column=1, sticky=Tk.E)

      # bol tail description
      def btailcontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         source=self.msource.get()
         if len(source) == 0 or source == 'None': source=None
         self.ztfp.data[objid].show_corner(
            ax, engine='bol_tail', model=self.tailmodel.get(), source=source,
            clobber=clobber,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.bbtop1, text="Contour", command=btailcontours,fg='blue')
      inputFrameBut.grid(row=3, column=1, sticky=Tk.SE)

      # bol full description
      def bfulldescription():
         model=self.bfmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.bbtop1, text="describe", command=bfulldescription, fg='brown')
      inputFrameBut.grid(row=4, column=1, sticky=Tk.E)
      
      # bol full contour
      def bfullcontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         source=self.msource.get()
         if len(source) == 0 or source == 'None': source=None
         self.ztfp.data[objid].show_corner(
            ax, engine='bol_full', model=self.bfmodel.get(), source=source,
            clobber=clobber,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.bbtop1, text="Contour", command=bfullcontours,fg='brown')
      inputFrameBut.grid(row=4, column=1, sticky=Tk.SE)

      # set texp
      Frame = Tk.LabelFrame(self.bbtop1, text="set t")
      Frame.grid(row=5,column=0,sticky=Tk.W)
      
      def sett0(): self.ztfp.data[objid].set_peak_bol_main()
      inputFrameBut = Tk.Button(Frame, text="t0", command=sett0)
      inputFrameBut.grid(row=0, column=0, sticky=Tk.NW)
      
      def settexp():
         self.ztfp.data[objid].set_texp_bol_main()
      inputFrameBut = Tk.Button(Frame, text="texp", command=settexp)
      inputFrameBut.grid(row=0, column=1, sticky=Tk.NW) 
      
   def spectrafit_popup(self, objid):
      # pop up window
      self.spectop = Tk.Toplevel(self.parent)
      self.spectop.geometry("500x500+50+50")      
      self.spectop.title('spectra fitter for %s'%objid)
      self.spectop.focus_set()
      
      ### specify model
      menu = Tk.LabelFrame(self.spectop, text="Spectral line fitter")           
      menu.grid_rowconfigure(0, weight=1)
      menu.grid_columnconfigure(0, weight=1)
      menu.grid_columnconfigure(1, weight=1)
      menu.grid_columnconfigure(2, weight=100)
      menu.grid(row=0, column=0)
      
      # data infos
      row1 = Tk.Frame(self.spectop)      
      row1.grid(row=0,column=0,sticky=Tk.W)      
      labelFrame = Tk.LabelFrame(row1, text="Spectral Fits")
      labelFrame.grid(row=0,column=0)
      
      self.datainfo5 = Tk.StringVar()
      self.datainfo5.set(self.checkdata(objid, ['sed']))
      tb = Tk.Label(labelFrame,textvariable=self.datainfo5, fg="red", justify= Tk.LEFT)
      tb.grid(row=0,column=0)         
      
      # MCMC
      # row1 right
      Frame = Tk.Frame(self.spectop)
      Frame.grid(row=0,column=1,sticky=Tk.W)

      # sources            
      _Frame = Tk.LabelFrame(Frame, text="Source")
      _Frame.grid(row=0,column=0,sticky=Tk.W)
      sources = [None]
      if 'spec' in self.ztfp.data[objid].__dict__:
         sources = np.append(sources, self.ztfp.data[objid].spec.sort_spectra(returnepochs=True))         
      self.specsource = Tk.StringVar()
      self.specsource.set( sources[0] )
      
      lcsourceMenu = Tk.OptionMenu(_Frame, self.specsource, *sources)
      lcsourceMenu.config(justify=Tk.LEFT,width=4)
      lcsourceMenu.grid(row=0,column=0,sticky=Tk.E)
      
      def updatesource(k): self.specsource.set( k )
      def updateinfos():
         sources = [None]
         if 'spec' in self.ztfp.data[objid].__dict__:
            sources = np.append(sources, list(self.ztfp.data[objid].spec.data.keys()))
         self.specsource.set( sources[0] )         
         # update option menu
         menu = lcsourceMenu["menu"]
         menu.delete(0, "end")
         for t in sources:            
            menu.add_command(label=t,command=lambda k=t: updatesource(k))  
      button = Tk.Button(Frame,text="update",command=updateinfos)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=1)
      
      def deletefit():
         specsource = self.specsource.get()
         specline = self.specline.get()
         if specline=='sntype':
            specline = constants.line_forsn[self.ztfp.data[objid].sntype]
         source = '%s_%s' % ('_'.join(specsource.split()), specline)
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'specline' in self.ztfp.data[objid].__dict__['fitcls']:
               try: del self.ztfp.data[objid].__dict__['fitcls']['specline'][source]
               except: pass
         self.datainfo5.set(self.checkdata(objid, ['sed']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo5, fg="red", justify= Tk.LEFT)         
      button = Tk.Button(Frame,text="del",command=deletefit)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=1,column=1)
      
      def reset():
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'specline' in self.ztfp.data[objid].__dict__['fitcls']:
               del self.ztfp.data[objid].__dict__['fitcls']['specline']
            if 'specv_evolution' in self.ztfp.data[objid].__dict__['fitcls']:
               del self.ztfp.data[objid].__dict__['fitcls']['specv_evolution']
         self.datainfo5.set(self.checkdata(objid, ['sed']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo5, fg="red", justify= Tk.LEFT)         
      button = Tk.Button(Frame,text="reset",command=reset)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=2,column=1)

      _Frame = Tk.LabelFrame(Frame, text="element")
      _Frame.grid(row=1,column=0)
      lines = ['full', 'sntype'] + list(constants.line_location.keys())
      self.specline = Tk.StringVar()
      self.specline.set( lines[0] )
      
      Menu = Tk.OptionMenu(_Frame, self.specline, *lines)
      Menu.config(justify=Tk.LEFT,width=4)
      Menu.grid(row=0,column=0,sticky=Tk.W)
      
      # spectra plots            
      Frame = Tk.LabelFrame(self.spectop, text="plotter")
      Frame.grid(row=2,column=0,columnspan=2,sticky=Tk.W)

      _Frame = Tk.LabelFrame(Frame, text="Type")
      _Frame.grid(row=0,column=0)
      _list = ['original', 'rest', 'bin', 'continuum', 'flat']
      self.spectype = Tk.StringVar()
      self.spectype.set( 'flat' )
      Menu = Tk.OptionMenu(_Frame, self.spectype, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      pfactor = snobject_pars['pfactor']
      l = Tk.Label(Frame,text='pfactor=',fg="black")      
      l.grid(row=0,column=1)
      self.pfactentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.pfactentry.insert(0, pfactor)
      self.pfactentry.grid(row=0, column=2)
      
      _Frame = Tk.LabelFrame(Frame, text="bin method")
      _Frame.grid(row=1,column=0)
      _list = ['sum', 'average', 'gauss', 'median', 'savgol']
      self.specbinmethod = Tk.StringVar()
      self.specbinmethod.set( 'savgol' )
      Menu = Tk.OptionMenu(_Frame, self.specbinmethod, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0, column=0, sticky=Tk.W)

      binsize = snobject_pars['bin_size']
      l = Tk.Label(Frame,text='bin_size=',fg="black")      
      l.grid(row=1,column=1)
      self.specbinsizeentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.specbinsizeentry.insert(0, binsize)
      self.specbinsizeentry.grid(row=1, column=2)
      
      _Frame = Tk.LabelFrame(Frame, text="continuum method")
      _Frame.grid(row=2,column=0)
      _list = ["scalar", "linear", "quadratic", "cubic", "poly", "exponential"]
      self.speccontmethod = Tk.StringVar()
      self.speccontmethod.set( 'poly' )
      Menu = Tk.OptionMenu(_Frame, self.speccontmethod, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      cdeg = snobject_pars['continuum_degree']
      l = Tk.Label(Frame,text='poly degree=',fg="black")      
      l.grid(row=2,column=1)
      self.speccontsizeentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.speccontsizeentry.insert(0, cdeg)
      self.speccontsizeentry.grid(row=2, column=2)
      
      def showspec():         
         if not 'spec' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse spectra for %s first'%objid)
            return                              
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         ax = self.ztfp.data[objid].ax1 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         
         source = self.specsource.get()
         if source == 'None': source=eval(source)
         else: source = [source]
         stype = self.spectype.get()
         if stype == 'None': stype=eval(stype)
         pfactor = self.pfactentry.get()
         bin_method = self.specbinmethod.get()
         bin_size = self.specbinsizeentry.get()
         continuum_method = self.speccontmethod.get()
         continuum_degree = self.speccontsizeentry.get()
         sn_line = self.specline.get()
         if sn_line == 'None': sn_line=eval(sn_line)
         
         # update spectra
         for _ in self.ztfp.data[objid].spec.data:            
            self.ztfp.data[objid].spec.data[_]['data']._bin_spectrum(
               bin_method=bin_method, bin_size=int(bin_size),
               savgol_order=int(snobject_pars['savgol_order']),
            )
            self.ztfp.data[objid].spec.data[_]['data']._estimate_continuum(
               continuum_method=continuum_method, continuum_degree=int(continuum_degree),
            )
            self.ztfp.data[objid].spec.data[_]['data']._correct_continuum() 
         # plot
         self.currentfigure = 'spectra'         
         self.ztfp.data[objid]._ax1(
            plot_specsources=source, show_stype=stype, show_element=sn_line,
            pfactor=float(pfactor), continuum_method=continuum_method,
            continuum_degree=int(continuum_degree),
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show spectra for %s'%objid)
         
      inputFrameBut = Tk.Button(self.spectop, text="Show", command=showspec)
      inputFrameBut.grid(row=2, column=2, sticky=Tk.E)
      
      # spectra fits
      Frame = Tk.LabelFrame(self.spectop, text="Line fit",fg="green")
      Frame.grid(row=3,column=0,columnspan=2,sticky=Tk.W)      
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=0)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.specroutine = Tk.StringVar()
      self.specroutine.set( snobject_pars['bol_early_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.specroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=0)      
      _list = get_model(which_engine='specline', with_alias=False)['specline']
      self.specmodel = Tk.StringVar()
      self.specmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.specmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      _Frame = Tk.LabelFrame(Frame, text="Forced range")
      _Frame.grid(row=2,column=0)      
      l = Tk.Label(_Frame,text='-',fg="black")      
      l.grid(row=0,column=1)
      self.lineminloc = Tk.Entry(_Frame, justify=Tk.LEFT, width=3)
      self.linemaxloc = Tk.Entry(_Frame, justify=Tk.RIGHT, width=3)
      self.lineminloc.insert(0, 'None')
      self.linemaxloc.insert(0, 'None')
      self.lineminloc.grid(row=0, column=0)
      self.linemaxloc.grid(row=0, column=2)
      
      #
      nsteps = snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=0,column=1)
      self.nstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepsentry.insert(0, nsteps)
      self.nstepsentry.grid(row=0, column=2)
      
      nsteps_burnin = snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=1,column=1)
      self.nstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepburnsentry.insert(0, nsteps_burnin)
      self.nstepburnsentry.grid(row=1, column=2)

      nwalkers = snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=2,column=1)
      self.nwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nwalkerentry.insert(0, nwalkers)
      self.nwalkerentry.grid(row=2, column=2)
      
      # done
      def runfit():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.specsource.get()
         sn_line = self.specline.get()
         if sn_line == 'None':
            self.Message('Define a line first')
            return
         self.Message('%s line fits for source: %s with clobber=%s'%(sn_line,source,clobber))
         if len(source) == 0 or source == 'None': source=None
         if len(sn_line) == 0 or sn_line == 'None': source=None

         if is_number(self.lineminloc.get()) and is_number(self.linemaxloc.get()):
            force_range = [float(self.lineminloc.get()), float(self.linemaxloc.get())]
         else:
            force_range = None
         self.ztfp.data[objid].run_fit(
            'specline',
            fit_methods = [self.specmodel.get()],            
            specline_routine = self.specroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            fit_redo = clobber,
            source = source,
            sn_line = sn_line,
            force_range = force_range,
         )
         self.datainfo5.set(self.checkdata(objid, ['sed']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo5, fg="red", justify= Tk.LEFT)         
      button2 = Tk.Button(self.spectop,text="Run",command=runfit,fg="green")
      button2.grid(row=3, column=2, sticky=Tk.NE)

      # description
      def specdescription():
         model=self.specmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.spectop, text="describe", command=specdescription, fg='green')
      inputFrameBut.grid(row=3, column=2, sticky=Tk.E)
      
      # corner plots
      def speccontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         _source=self.specsource.get()          
         _sline = self.specline.get()         
         if _sline=='full':
            self.Message('Warning: specific a line first')
            return            
         elif _sline=='sntype': _sline = constants.line_forsn[self.ztfp.data[objid].sntype]        
         else: pass 
         if _source != 'None' and _sline != 'None': 
            source = '%s_%s' % ('_'.join(_source.split()), _sline)
         else:
            source=None
         self.ztfp.data[objid].show_corner(
            ax, engine='specline', model=self.specmodel.get(), source=source,
            clobber=clobber,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)
      inputFrameBut = Tk.Button(self.spectop, text="Contours", command=speccontours,fg="green")
      inputFrameBut.grid(row=3, column=2, sticky=Tk.SE)
      
      # line velocity      
      Frame = Tk.LabelFrame(self.spectop, text="Line velocity fit",fg="blue")
      Frame.grid(row=4,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['specv_evolution_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")
      l.grid(row=0,column=0)
      self.specvminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specvminentry.insert(0, gpmin)
      self.specvminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.specvmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specvmaxentry.insert(0, gpmax)
      self.specvmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['specv_evolution_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.specvpminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specvpminentry.insert(0, gpmin)
      self.specvpminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.specvpmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specvpmaxentry.insert(0, gpmax)
      self.specvpmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.specvroutine = Tk.StringVar()
      self.specvroutine.set( snobject_pars['specv_evolution_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.specvroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)      
      _list = get_model(which_engine='specv_evolution', with_alias=False)['specv_evolution']
      self.specvmodel = Tk.StringVar()
      self.specvmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.specvmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)      

      def runlv():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         fit_fitr = [float(self.specvminentry.get()),float(self.specvmaxentry.get())]
         fit_fitp = [float(self.specvpminentry.get()),float(self.specvpmaxentry.get())]         
         self.ztfp.data[objid].run_fit(
            'specv_evolution',
            fit_methods = [self.specvmodel.get()],
            specv_evolution_xrange = fit_fitr,
            specv_evolution_xrangep = fit_fitp,
            specv_evolution_routine = self.specvroutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),            
            fit_redo=clobber,
         )
         self.datainfo5.set(self.checkdata(objid, ['sed']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo5, fg="red", justify= Tk.LEFT)      
      button1 = Tk.Button(self.spectop,text="Run",command=runlv,fg='blue')
      button1.grid(row=4, column=2, sticky=Tk.NE)

      def showlv():
         sn_line = self.specline.get()
         if sn_line == 'None':
            self.Message('Define a line first')
            return         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'baseline'
         ax = self.ztfp.data[objid].ax6 = self.add_subplot(111, zoomable='horizontal', cursor='both')         
         self.ztfp.data[objid]._ax6(
            sn_line = sn_line,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Spectral velocity evolution plot for %s'%objid)    
      button1 = Tk.Button(self.spectop,text="Show",command=showlv,fg='blue')
      button1.grid(row=5, column=2, sticky=Tk.NE)
      
      # description
      def lvdescription():
         model=self.specvmodel.get()
         infos = get_pars(model, with_alias=True)
         for n, k in enumerate(['enginename', 'parname', 'par', 'same_parameter',
                                'fit_error', 'bestv', 'bounds']):
            if n == 0:
               self.Message('model: %s\n%s: %s\n' % (model, k, infos[k]),disable=False,restart=True)
            else:
               self.Message('%s: %s\n' % (k, infos[k]),disable=False,restart=False)            
      inputFrameBut = Tk.Button(self.spectop, text="describe", command=lvdescription, fg='blue')
      inputFrameBut.grid(row=4, column=2, sticky=Tk.E)
      
      # bol early corner
      def lvcontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both') 
         self.ztfp.data[objid].show_corner(
            ax, engine='specv_evolution', model=self.specvmodel.get(), source=None,
            clobber=clobber,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.spectop, text="Contour", command=lvcontours,fg='blue')
      inputFrameBut.grid(row=4, column=2, sticky=Tk.SE)
      
   def baselinetop_popup(self, objid):
      #objid = self.obj.get()
      
      # pop up window
      self.baselinetop = Tk.Toplevel(self.parent)
      self.baselinetop.geometry("200x200+50+50")      
      self.baselinetop.title('Check baseline')

      xmin = -100
      l = Tk.Label(self.baselinetop,text='xmin=',fg="black")      
      l.grid(row=0,column=0)
      self.xminentry = Tk.Entry(self.baselinetop, justify=Tk.LEFT, width=6)
      self.xminentry.insert(0, xmin)
      self.xminentry.grid(row=0, column=1)

      xmax = -20
      l = Tk.Label(self.baselinetop,text='xmax=',fg="black")      
      l.grid(row=1,column=0)
      self.xmaxentry = Tk.Entry(self.baselinetop, justify=Tk.LEFT, width=6)
      self.xmaxentry.insert(0, xmax)
      self.xmaxentry.grid(row=1, column=1)      
      
      def showbaseline():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'baseline'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')         
         baseline = self.ztfp.data[objid].calibrate_baseline(ax=ax, key='fcqfid', source='ztffp',
               xmin=float(self.xminentry.get()), xmax=float(self.xmaxentry.get()),
               ax_xlim=None, ax_ylim=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Baseline plots for %s'%objid)
      inputFrameBut = Tk.Button(self.baselinetop, text="Show", command=showbaseline)
      inputFrameBut.grid(row=2, column=0, sticky=Tk.W)

      def correctbaseline():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'baseline'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         baseline = self.ztfp.data[objid].calibrate_baseline(ax=None, key='fcqfid', source='ztffp',
               xmin=float(self.xminentry.get()), xmax=float(self.xmaxentry.get()),
               ax_xlim=None, ax_ylim=None)
         self.ztfp.data[objid].correct_baseline(baseline, key='fcqfid', source='ztffp')
         baseline = self.ztfp.data[objid].calibrate_baseline(ax=ax, key='fcqfid', source='ztffp',
               xmin=float(self.xminentry.get()), xmax=float(self.xmaxentry.get()),
               ax_xlim=None, ax_ylim=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Baseline plots for %s'%objid)
      inputFrameBut = Tk.Button(self.baselinetop, text="Correct", command=correctbaseline)
      inputFrameBut.grid(row=3, column=0, sticky=Tk.W)

   ### popup for sed plot ###   
   def sed_popup(self, objid):
      
      # pop up window
      self.sedtop = Tk.Toplevel(self.parent)
      self.sedtop.geometry("500x500+50+50")      
      self.sedtop.title('Spectral Energy Distribution for %s'%objid)
      
      # row1 
      row1 = Tk.Frame(self.sedtop)
      row1.grid(row=0,column=0,sticky=Tk.W)      
      labelFrame = Tk.LabelFrame(row1, text="Fits")
      labelFrame.grid(row=0,column=0)
      
      self.seddatainfo = Tk.StringVar()
      self.seddatainfo.set(self.checkdata(objid, ['sedfits']))
      tb = Tk.Label(labelFrame,textvariable=self.seddatainfo, fg="red", justify= Tk.LEFT)
      tb.grid(row=0,column=0)
      
      # MCMC
      # row1 right
      Frame = Tk.Frame(self.sedtop)
      Frame.grid(row=0,column=1,sticky=Tk.W)
      
      nsteps = snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=0,column=0)
      self.sednstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.sednstepsentry.insert(0, nsteps)
      self.sednstepsentry.grid(row=0, column=1)
      
      nsteps_burnin = snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=1,column=0)
      self.sednstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.sednstepburnsentry.insert(0, nsteps_burnin)
      self.sednstepburnsentry.grid(row=1, column=1)

      nwalkers = snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=0,column=2)
      self.sednwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.sednwalkerentry.insert(0, nwalkers)
      self.sednwalkerentry.grid(row=0, column=3)
            
      phase = 0
      l = Tk.Label(Frame,text='phase=',fg="black")      
      l.grid(row=1,column=2)
      self.sedphaseentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.sedphaseentry.insert(0, phase)
      self.sedphaseentry.grid(row=1, column=3)
      
      def reset():         
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'sed' in self.ztfp.data[objid].__dict__['fitcls']:
               del self.ztfp.data[objid].__dict__['fitcls']['sed']            
         self.seddatainfo.set(self.checkdata(objid, ['sedfits']))
         tb = Tk.Label(labelFrame,textvariable=self.seddatainfo, fg="red", justify= Tk.LEFT)
      button = Tk.Button(Frame,text="reset",command=reset)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=2,column=2, sticky=Tk.SW)
      
      # multiband BB fits         
      Frame = Tk.LabelFrame(self.sedtop, text="Multiband BB fit",fg="blue")
      Frame.grid(row=1,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['sed_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.sedfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sedfitminentry.insert(0, gpmin)
      self.sedfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.sedfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sedfitmaxentry.insert(0, gpmax)
      self.sedfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['sed_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.sedplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sedplotminentry.insert(0, gpmin)
      self.sedplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.sedplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.sedplotmaxentry.insert(0, gpmax)
      self.sedplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.sedroutine = Tk.StringVar()
      self.sedroutine.set( snobject_pars['sed_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.sedroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = get_model(which_engine='sed', with_alias=False)['sed']
      self.sedmodel = Tk.StringVar()
      self.sedmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.sedmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      if 'lc' in self.ztfp.data[objid].__dict__:
         colours = list(np.unique(self.ztfp.data[objid].lc['filter']))
      else:
         colours = snobject_pars['sed_bands']
      l = Tk.Label(Frame,text='colours=',fg="black")      
      l.grid(row=2,column=0)
      self.sedcolorentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.sedcolorentry.insert(0, colours)
      self.sedcolorentry.grid(row=2, column=1)

      self.sedbbmixfit = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(Frame, text="fit together",variable=self.sedbbmixfit) 
      cButton.grid(row=2,column=3,columnspan=2,sticky=Tk.W)
      
      # spectral BB fits         
      Frame = Tk.LabelFrame(self.sedtop, text="spectra BB fit",fg="green")
      Frame.grid(row=2,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['sed_xrange']      
      l = Tk.Label(Frame,text='Fit range:',fg="black")      
      l.grid(row=0,column=0)
      self.specsedfitminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specsedfitminentry.insert(0, gpmin)
      self.specsedfitminentry.grid(row=0, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=0,column=1, sticky=Tk.E)
      self.specsedfitmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specsedfitmaxentry.insert(0, gpmax)
      self.specsedfitmaxentry.grid(row=0, column=2, sticky=Tk.W)
      
      gpmin, gpmax = snobject_pars['sed_xrangep']      
      l = Tk.Label(Frame,text='Plot range:',fg="black")      
      l.grid(row=1,column=0)
      self.specsedplotminentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specsedplotminentry.insert(0, gpmin)
      self.specsedplotminentry.grid(row=1, column=1, sticky=Tk.W)
      l = Tk.Label(Frame,text='-',fg="black")      
      l.grid(row=1,column=1, sticky=Tk.E)
      self.specsedplotmaxentry = Tk.Entry(Frame, justify=Tk.LEFT, width=4)
      self.specsedplotmaxentry.insert(0, gpmax)
      self.specsedplotmaxentry.grid(row=1, column=2, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=3)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.specsedroutine = Tk.StringVar()
      self.specsedroutine.set( snobject_pars['sed_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.specsedroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = get_model(which_engine='sed', with_alias=False)['sed']
      self.specsedmodel = Tk.StringVar()
      self.specsedmodel.set( _list[0] )
      Menu = Tk.OptionMenu(_Frame, self.specsedmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      assert 'lc' in self.ztfp.data[objid].__dict__      
      speccolours = list(np.unique(self.ztfp.data[objid].lc['filter']))
      speccolour = speccolours[0]

      l = Tk.Label(Frame,text='abs cal band=',fg="black")      
      l.grid(row=2,column=0)
      self.sedspeccolor = Tk.StringVar()
      self.sedspeccolor.set( speccolour )
      Menu = Tk.OptionMenu(Frame, self.sedspeccolor, *speccolours)
      Menu.config(justify=Tk.LEFT,width=2)
      Menu.grid(row=2,column=1, sticky=Tk.W)

      l = Tk.Label(Frame,text='abs cal method=',fg="black")      
      l.grid(row=2,column=3)
      self.sedspecmethod = Tk.StringVar()
      self.sedspecmethod.set( 'cal' )
      Menu = Tk.OptionMenu(Frame, self.sedspecmethod, *['cal', 'mangle'])
      Menu.config(justify=Tk.LEFT,width=2)
      Menu.grid(row=2,column=4, sticky=Tk.W)   
      
      def makesed1():
         bbphase = float(self.sedphaseentry.get())
         cinterp = []
         for index in self.interplist.curselection():
            cinterp.append(self.interplist.get(index))
         if len(cinterp) == 0:
            self.Message('Select at least one option to match colours')
            return
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         fit_bands = self.sedcolorentry.get().split()
         
         # calculate colours
         self.ztfp.data[objid].bb_colors(
            xpred=[bbphase],
            sed_bands=fit_bands,
            sed_color_interp=cinterp,
            tdbin = int(self.tdbentry1.get()),            
            snrt = int(self.snrentry1.get()),
            corr_mkw = corr_mkw,
            corr_host = corr_host,
            index=0,
            returnv=False
         )
         
         # run fitting
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         if self.sedbbmixfit.get()==0: sed_color_interp_mix=False
         else: sed_color_interp_mix=True
         fit_fitr = [float(self.sedfitminentry.get()),float(self.sedfitmaxentry.get())]
         fit_fitp = [float(self.sedplotminentry.get()),float(self.sedplotmaxentry.get())]   
         self.ztfp.data[objid].run_fit(
            'sed',
            source='bb',
            make_bol=['bb'],
            fit_methods = [self.sedmodel.get()],
            sed_bands = fit_bands,
            sed_color_interp=cinterp,
            sed_xrange = fit_fitr,
            sed_xrangep = fit_fitp,
            sed_routine = self.sedroutine.get(),
            nwalkers = int(self.sednwalkerentry.get()),
            nsteps = int(self.sednstepsentry.get()),
            nsteps_burnin = int(self.sednstepburnsentry.get()),
            fit_redo = clobber,
            sed_color_interp_mix = sed_color_interp_mix,
         )
         
         # calculate bolometric LC
         self.ztfp.data[objid].bb_bol(
            make_bol = ['bb'],
            fastsedfitting=False,
            sed_bands = fit_bands,
            sed_color_interp=cinterp,
            sed_xrangep = fit_fitp,
         )
         
         # update info
         self.seddatainfo.set(self.checkdata(objid, ['sedfits']))
         tb = Tk.Label(labelFrame,textvariable=self.seddatainfo, fg="red", justify= Tk.LEFT)
         
      data1button = Tk.Button(self.sedtop,text="run",command=makesed1)
      data1button.grid(row=1, column=2, sticky=Tk.NE)

      def plotsed1():
         bbphase = float(self.sedphaseentry.get())
         cinterp = []
         for index in self.interplist.curselection():
            cinterp.append(self.interplist.get(index))
         if len(cinterp) == 0: cinterp = None
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         fit_bands = self.sedcolorentry.get().split()
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'BB'                             
         ax = self.ztfp.data[objid].ax5 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax5(
            bbphase, index=0,
            interpolation = None,
            corr_mkw = corr_mkw,
            corr_host = corr_host,
            tdbin = int(self.tdbentry1.get()),            
            snrt = int(self.snrentry1.get()),
            do_kcorr = True,
            ab2vega = True,
            sed_type = ['bb'],
            sed_color_interp = cinterp,
            sed_bands = fit_bands,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('BB plot')
      data1button = Tk.Button(self.sedtop,text="show",command=plotsed1)      
      data1button.grid(row=1, column=2, sticky=Tk.SE)

      def sedbbcontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         
         bbphase = float(self.sedphaseentry.get())
         if self.sedbbmixfit.get()==0: sed_color_interp_mix=False
         else: sed_color_interp_mix=True
         jd = self.ztfp.data[objid].t0 + bbphase*(1+self.ztfp.data[objid].z)
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass         
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         if sed_color_interp_mix:
            sourcename = 'bb_mix_%.1f' % (jd)
         else:
            interpolation = self.interplist.get(Tk.ACTIVE)            
            sourcename = 'bb_%s_%.1f' % (interpolation, jd)
         self.ztfp.data[objid].show_corner(
            ax, engine='sed', source=sourcename,
            model=self.sedmodel.get(), clobber=clobber, filts=None
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.sedtop, text="Contours", command=sedbbcontours)
      inputFrameBut.grid(row=1, column=2, sticky=Tk.E)
      
      def makesed2():
         if 'spec' not in self.ztfp.data[objid].__dict__:
            self.Message('Warning: parse spctra first')
            return
         bbphase = float(self.sedphaseentry.get())         
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         
         # run fitting
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True         
         fit_fitr = [float(self.sedfitminentry.get()),float(self.sedfitmaxentry.get())]
         fit_fitp = [float(self.sedplotminentry.get()),float(self.sedplotmaxentry.get())]   
         self.ztfp.data[objid].run_fit(
            'sed',
            source='spec',
            make_bol=['spec'],
            fit_methods = [self.sedmodel.get()],            
            sed_xrange = fit_fitr,
            sed_xrangep = fit_fitp,
            sed_routine = self.sedroutine.get(),
            nwalkers = int(self.sednwalkerentry.get()),
            nsteps = int(self.sednstepsentry.get()),
            nsteps_burnin = int(self.sednstepburnsentry.get()),
            fit_redo = clobber,
            sed_abscal_bands = [self.sedspeccolor.get()],
            sed_abscal_method = self.sedspecmethod.get(),
            specfit_phase = [bbphase-int(self.tdbentry1.get()),
                             bbphase+int(self.tdbentry1.get())],
         )
         
         # calculate bolometric LC
         self.ztfp.data[objid].bb_bol(
            make_bol = ['spec'],
            fastsedfitting=False,
            sed_xrangep = fit_fitp,
         )
         
         # update info
         self.seddatainfo.set(self.checkdata(objid, ['sedfits']))
         tb = Tk.Label(labelFrame,textvariable=self.seddatainfo, fg="red", justify= Tk.LEFT)
      data1button = Tk.Button(self.sedtop,text="run",command=makesed2)
      data1button.grid(row=2, column=2, sticky=Tk.NE)

      def sedspeccontours():
         #if self.clobberVar.get()==0: clobber=False
         #else: clobber=True
         clobber = True
         bbphase = float(self.sedphaseentry.get())
         
         pspec = None
         for _ in self.ztfp.data[objid].spec.data:           
            _phase  = float(self.ztfp.data[objid].spec.data[_]['phase'])
            if abs(_phase - bbphase) > int(self.tdbentry1.get()): continue
            pspec = _phase
         if pspec is None:
            self.Message('No spectra found at %.2f within %i day bin' % (bbphase, int(self.tdbentry1.get())))
            return
         jd = pspec * (1+self.ztfp.data[objid].z) + self.ztfp.data[objid].t0
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass         
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')                 
         sourcename = 'spec_%.1f' % (jd)
         self.ztfp.data[objid].show_corner(
            ax, engine='sed', source= sourcename,
            model=self.sedmodel.get(), clobber=clobber, filts=None
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.sedtop, text="Contours", command=sedspeccontours)
      inputFrameBut.grid(row=2, column=2, sticky=Tk.E)
      
      def plotsed2():
         bbphase = float(self.sedphaseentry.get())         
         if self.corrmkwebv.get() == 0: corr_mkw=False
         else:  corr_mkw=True
         if self.corrhostebv.get() == 0: corr_host=False
         else:  corr_host=True
         
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'BB'                             
         ax = self.ztfp.data[objid].ax5 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax5(
            bbphase, 
            corr_mkw = corr_mkw,
            corr_host = corr_host,
            tdbin = int(self.tdbentry1.get()),            
            snrt = int(self.snrentry1.get()),
            do_kcorr = True,
            ab2vega = True,
            sed_type = ['spec'],
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('BB plot')
      data1button = Tk.Button(self.sedtop,text="show",command=plotsed2)      
      data1button.grid(row=2, column=2, sticky=Tk.SE)
      
   ### popup for Flux LC ###  
   def fluxfit_popup(self, objid):
      
      # pop up window
      self.fluxtop = Tk.Toplevel(self.parent)
      self.fluxtop.geometry("500x500+50+50")      
      self.fluxtop.title('Flux model generator for %s'%objid)

      # sources            
      _Frame = Tk.LabelFrame(self.fluxtop, text="Source")
      _Frame.grid(row=0,column=0,sticky=Tk.W)
      sources = [None]
      if 'lc' in self.ztfp.data[objid].__dict__:
         sources = np.append(sources, np.unique(self.ztfp.data[objid].lc['source']))
      self.lcsource1 = Tk.StringVar()
      self.lcsource1.set( sources[0] )
      
      lcsourceMenu = Tk.OptionMenu(_Frame, self.lcsource1, *sources)
      lcsourceMenu.config(justify=Tk.LEFT,width=4)
      lcsourceMenu.grid(row=0,column=0,sticky=Tk.E)

      # filters            
      _Frame = Tk.LabelFrame(self.fluxtop, text="Filter")
      _Frame.grid(row=1,column=0,sticky=Tk.W)
      filters = [None]
      if 'lc' in self.ztfp.data[objid].__dict__:
         filters = np.append(filters, np.unique(self.ztfp.data[objid].lc['filter']))      
      self.lcfs = Tk.StringVar()
      self.lcfs.set( filters[0] )
      
      lcfMenu = Tk.OptionMenu(_Frame, self.lcfs, *filters)
      lcfMenu.config(justify=Tk.LEFT,width=4)
      lcfMenu.grid(row=0,column=0,sticky=Tk.E)

      # x style
      _Frame = Tk.LabelFrame(self.fluxtop, text="x style")
      _Frame.grid(row=0,column=0)
      self.lcxstyle = Tk.StringVar()
      self.lcxstyle.set( 'jd' )
      lcxMenu = Tk.OptionMenu(_Frame, self.lcxstyle, *['jd', 'rp'])
      lcxMenu.grid(row=0,column=0)
      
      # y style
      _Frame = Tk.LabelFrame(self.fluxtop, text="y style")
      _Frame.grid(row=1,column=0)
      self.lcystyle = Tk.StringVar()
      self.lcystyle.set( 'original' )
      lcyMenu = Tk.OptionMenu(_Frame, self.lcystyle, *['original', 'norm'])
      lcyMenu.grid(row=0,column=0)
      
      # show GP
      self.showgpVar = Tk.IntVar(value=0) 
      showgpButton = Tk.Checkbutton(self.fluxtop, text="Show GP",variable=self.showgpVar)      
      showgpButton.config(justify=Tk.LEFT,width=12)
      showgpButton.grid(row=0,column=0,sticky=Tk.E)      
      
      # show fits
      self.showfitsVar = Tk.IntVar(value=0) 
      showfitsutton = Tk.Checkbutton(self.fluxtop, text="Show fits",variable=self.showfitsVar)      
      showfitsutton.config(justify=Tk.LEFT,width=12)
      showfitsutton.grid(row=1,column=0,sticky=Tk.E)      
      
      def updatesource(k): self.lcsource1.set( k )
      def updatefilter(k): self.lcfs.set( k )
      def updateinfos():
         sources = [None]
         if 'lc' in self.ztfp.data[objid].__dict__:            
            sources = np.append(sources, np.unique(self.ztfp.data[objid].lc['source']))
         self.lcsource1.set( sources[0] )

         filters = [None]
         if 'lc' in self.ztfp.data[objid].__dict__:
            filters = np.append(filters, np.unique(self.ztfp.data[objid].lc['filter']))
         self.lcfs.set( filters[0] )
         
         # update option menu
         menu = lcsourceMenu["menu"]
         menu.delete(0, "end")
         for t in sources:            
            menu.add_command(label=t,command=lambda k=t: updatesource(k))

         menu = lcfMenu["menu"]
         menu.delete(0, "end")
         for t in filters:
            menu.add_command(label=t,command=lambda k=t: updatefilter(k))
      button = Tk.Button(self.fluxtop,text="update",command=updateinfos)
      button.config(justify=Tk.LEFT,width=2)
      button.grid(row=2,column=0,sticky=Tk.W)
      
      def updateplot():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         source = self.lcsource1.get()
         filters = self.lcfs.get()
         if filters == 'None': filters = eval(filters)         
         if source == 'None': source = eval(source)
         showgp = self.showgpVar.get()
         if showgp == 0: show_gp=False
         else: show_gp=True
         showfits = self.showfitsVar.get()
         if showfits == 0: show_fit=False
         else: show_fit=True                  
         self.currentfigure = 'Flux LC'
         ax = self.ztfp.data[objid].ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax(
            show_title=False, plot_sources=source, plot_bands=filters,
            show_fits=show_fit, show_gp=show_gp, ax_xstyle=self.lcxstyle.get(),
            ax_ystyle=self.lcystyle.get(),
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show Flux lcs for %s'%objid)
         
      button = Tk.Button(self.fluxtop,text="show",command=updateplot)
      button.config(justify=Tk.LEFT,width=2)
      button.grid(row=2,column=1, sticky=Tk.E)
      
      # Bazin fits
      menu = Tk.LabelFrame(self.fluxtop, text="Bazin model")      
      menu.grid(row=3, column=0) 
      
      # default p
      filt = 'r'
      a_p = 100
      dt_p = 0
      trise_p = 10
      tfall_p = 10
      c_p = 60            
      if objid is not None:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'multiband_main' in self.ztfp.data[objid].fitcls:
               for model in self.ztfp.data[objid].fitcls['multiband_main']:
                  if model not in ['bazin1', 'bazin2']: continue
                  try:a_p = self.ztfp.data[objid].fitcls['multiband_main'][model].get_par(filt=filt, parname='a')
                  except:pass
                  try:trise_p = self.ztfp.data[objid].fitcls['multiband_main'][model].get_par(filt=filt, parname='trise')
                  except:pass
                  try:tfall_p = self.ztfp.data[objid].fitcls['multiband_main'][model].get_par(filt=filt, parname='tfall')
                  except:pass
                  try:c_p = self.ztfp.data[objid].fitcls['multiband_main'][model].get_par(filt=filt, parname='c')
                  except:pass
      # a
      l = Tk.Label(menu,text='A =',fg="black")      
      l.grid(row=0,column=0)                  
      self.var_a = Tk.DoubleVar(value=a_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_a, wrap=True,
                        width=10, command=self.show_bazinmodel)
      def update_a(value):
         self.var_a.set( float(value) )      
         self.show_bazinmodel()
      slide = Tk.Scale(menu, variable=self.var_a, orient='horizontal',
                       length=200, command=update_a)
      spin['to'] = 1000
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 1000
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=0, column=1, sticky='news')
      slide.grid(row=0, column=2, sticky='news')

      # dt
      l = Tk.Label(menu,text='t0 =',fg="black")      
      l.grid(row=1,column=0)                  
      self.var_dt = Tk.DoubleVar(value=dt_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_dt, wrap=True,
                        width=10, command=self.show_bazinmodel)
      def update_dt(value):
         self.var_dt.set( float(value) )      
         self.show_bazinmodel()
      slide = Tk.Scale(menu, variable=self.var_dt, orient='horizontal',
                       length=200, command=update_dt)
      spin['to'] = 100
      spin['from'] = -100
      spin['increment'] = 0.01
      slide['to'] = 100
      slide['from'] = -100
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=1, column=1, sticky='news')
      slide.grid(row=1, column=2, sticky='news')

      # taurise
      l = Tk.Label(menu,text='trise =',fg="black")      
      l.grid(row=2,column=0)
      self.var_trise = Tk.DoubleVar(value=trise_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_trise, wrap=True,
                        width=10, command=self.show_bazinmodel)
      def update_trise(value):
         self.var_trise.set( float(value) )      
         self.show_bazinmodel()
      slide = Tk.Scale(menu, variable=self.var_trise, orient='horizontal',
                       length=200, command=update_trise)
      spin['to'] = 60
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 60
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=2, column=1, sticky='news')
      slide.grid(row=2, column=2, sticky='news')
      
      # taufall
      l = Tk.Label(menu,text='tfall =',fg="black")      
      l.grid(row=3,column=0)                  
      self.var_tfall = Tk.DoubleVar(value=tfall_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_tfall, wrap=True,
                        width=10, command=self.show_bazinmodel)
      def update_tfall(value):
         self.var_tfall.set( float(value) )      
         self.show_bazinmodel()
      slide = Tk.Scale(menu, variable=self.var_tfall, orient='horizontal',
                       length=200, command=update_tfall)
      spin['to'] = 120
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 120
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=3, column=1, sticky='news')
      slide.grid(row=3, column=2, sticky='news')

      # c
      l = Tk.Label(menu,text='C =',fg="black")      
      l.grid(row=4,column=0)                  
      self.var_c = Tk.DoubleVar(value=c_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_c, wrap=True,
                        width=10, command=self.show_bazinmodel)
      def update_c(value):
         self.var_c.set( float(value) )      
         self.show_bazinmodel()
      slide = Tk.Scale(menu, variable=self.var_c, orient='horizontal',
                       length=200, command=update_c)
      spin['to'] = 1000
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 1000
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=4, column=1, sticky='news')
      slide.grid(row=4, column=2, sticky='news')      
      
      # show
      self.showbazinVar = Tk.IntVar(value=0) 
      showbazinButton = Tk.Checkbutton(menu, text="Show",variable=self.showbazinVar)      
      showbazinButton.grid(row=5,column=0,sticky=Tk.W)     

      # power law fits      
      menu = Tk.LabelFrame(self.fluxtop, text="power law model")      
      menu.grid(row=4, column=0)
      
      # default p
      filt = 'r'
      a_p = 100
      texp_p = -18
      alpha_p = 2      
      c_p = 60            
      if objid is not None:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'multiband_early' in self.ztfp.data[objid].fitcls:
               for model in self.ztfp.data[objid].fitcls['multiband_early']:
                  if model not in ['powerlaw_multiple', 'powerlaw_single']: continue                 
                  try:a_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='a')
                  except:pass
                  try:texp_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='texp')
                  except:pass
                  try:alpha_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='alpha')
                  except:pass
                  try:c_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='c')
                  except:pass
      # a
      l = Tk.Label(menu,text='A =',fg="black")      
      l.grid(row=0,column=0)                  
      self.var_a1 = Tk.DoubleVar(value=a_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_a1, wrap=True,
                        width=10, command=self.show_plmodel)
      def update_a1(value):
         self.var_a1.set( float(value) )      
         self.show_plmodel()
      slide = Tk.Scale(menu, variable=self.var_a1, orient='horizontal',
                       length=200, command=update_a1)
      spin['to'] = 200
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 200
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=0, column=1, sticky='news')
      slide.grid(row=0, column=2, sticky='news')

      # texp
      l = Tk.Label(menu,text='texp =',fg="black")      
      l.grid(row=1,column=0) 
      self.var_texp1 = Tk.DoubleVar(value=texp_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_texp1, wrap=True,
                        width=10, command=self.show_plmodel)
      def update_texp(value):
         self.var_texp1.set( float(value) )      
         self.show_plmodel()
      slide = Tk.Scale(menu, variable=self.var_texp1, orient='horizontal',
                       length=200, command=update_texp)
      spin['to'] = 0
      spin['from'] = -100
      spin['increment'] = 0.01
      slide['to'] = 0
      slide['from'] = -100
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=1, column=1, sticky='news')
      slide.grid(row=1, column=2, sticky='news')

      # alpha
      l = Tk.Label(menu,text='alpha =',fg="black")      
      l.grid(row=2,column=0) 
      self.var_alpha = Tk.DoubleVar(value=alpha_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_alpha, wrap=True,
                        width=10, command=self.show_plmodel)
      def update_alpha(value):
         self.var_alpha.set( float(value) )      
         self.show_plmodel()
      slide = Tk.Scale(menu, variable=self.var_alpha, orient='horizontal',
                       length=200, command=update_alpha)
      spin['to'] = 10
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 10
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=2, column=1, sticky='news')
      slide.grid(row=2, column=2, sticky='news')
      
      # c
      l = Tk.Label(menu,text='C =',fg="black")      
      l.grid(row=3,column=0)                  
      self.var_c1 = Tk.DoubleVar(value=c_p)
      spin = Tk.Spinbox(menu, textvariable=self.var_c1, wrap=True,
                        width=10, command=self.show_plmodel)
      def update_c1(value):
         self.var_c1.set( float(value) )      
         self.show_plmodel()
      slide = Tk.Scale(menu, variable=self.var_c1, orient='horizontal',
                       length=200, command=update_c1)
      spin['to'] = 100
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 100
      slide['from'] = 0
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=3, column=1, sticky='news')
      slide.grid(row=3, column=2, sticky='news')      
      
      # show
      self.showplVar = Tk.IntVar(value=0) 
      showplButton = Tk.Checkbutton(menu, text="show",variable=self.showplVar)      
      showplButton.grid(row=4,column=0,sticky=Tk.W)      
      
   ### popup for bolometric LC fitter ###   
   def mbolfit_popup(self):
      # pop up window
      self.mfittop = Tk.Toplevel(self.parent)
      self.mfittop.geometry("500x500+50+50")      
      self.mfittop.title('model mbol generator')
      #self.top.grab_set()
      #self.mfittop.focus_set()
      
      ### specify model
      menu = Tk.LabelFrame(self.mfittop, text="Bolometric LC model")           
      menu.grid_rowconfigure(0, weight=1)
      menu.grid_columnconfigure(0, weight=1)
      menu.grid_columnconfigure(1, weight=1)
      menu.grid_columnconfigure(2, weight=100)
      menu.grid(row=0, column=0)     
      
      # mni
      l = Tk.Label(menu,text='M Ni (Msun) =',fg="black")      
      l.grid(row=0,column=0)
      self.var_mni = Tk.DoubleVar(value=.2)
      spin = Tk.Spinbox(menu, textvariable=self.var_mni, wrap=True,
                        width=10, command=self.show_mbolmodel)
      def update_mni(value):
         self.var_mni.set( float(value) )      
         self.show_mbolmodel()
      slide = Tk.Scale(menu, variable=self.var_mni, orient='horizontal',
                       length=200, command=update_mni)
      spin['to'] = 5
      spin['from'] = .01
      spin['increment'] = 0.01
      slide['to'] = 5
      slide['from'] = .01
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=0, column=1, sticky='news')
      slide.grid(row=0, column=2, sticky='news')

      # mej
      def update_mej(value):
         self.var_mej.set( float(value) )
         self.var_taum.set( radioactivemodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()) )
         self.var_t0.set( radioactivemodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.var_vexp.set( radioactivemodels.Mej_Ek_to_vej(self.var_mej.get(), self.var_ek.get()) )     
         self.show_mbolmodel()
      l = Tk.Label(menu,text='M ej (Msun) =',fg="black")      
      l.grid(row=1,column=0)
      self.var_mej = Tk.DoubleVar(value=1.)
      spin = Tk.Spinbox(menu, textvariable=self.var_mej, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_mej, orient='horizontal',
                       length=200, command=update_mej)
      spin['to'] = 50.
      spin['from'] = .1
      spin['increment'] = 0.1
      slide['to'] = 50.
      slide['from'] = .1
      slide['digits'] = 3
      slide['resolution'] = 0.1      
      spin.grid(row=1, column=1, sticky='news')
      slide.grid(row=1, column=2, sticky='news')
      
      # ek
      def update_ek(value):
         self.var_ek.set( float(value) )      
         self.var_taum.set( radioactivemodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()) )
         self.var_t0.set( radioactivemodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.var_vexp.set( radioactivemodels.Mej_Ek_to_vej(self.var_mej.get(), self.var_ek.get()) )     
         self.show_mbolmodel()
      l = Tk.Label(menu,text='E kin (foe) =',fg="black")      
      l.grid(row=2,column=0)
      self.var_ek = Tk.DoubleVar(value=1.)
      spin = Tk.Spinbox(menu, textvariable=self.var_ek, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_ek, orient='horizontal',
                       length=200, command=update_ek)
      spin['to'] = 100.
      spin['from'] = .1
      spin['increment'] = 0.1
      slide['to'] = 100.
      slide['from'] = .1
      slide['digits'] = 3
      slide['resolution'] = 0.1
      spin.grid(row=2, column=1, sticky='news')
      slide.grid(row=2, column=2, sticky='news')
      
      # texp
      def update_texp(value):
         self.var_texp.set( float(value) )      
         self.show_mbolmodel()
      l = Tk.Label(menu,text='T exp (day) =',fg="black")      
      l.grid(row=3,column=0)
      self.var_texp = Tk.DoubleVar(value=-18)
      spin = Tk.Spinbox(menu, textvariable=self.var_texp, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_texp, orient='horizontal',
                       length=200, command=update_texp)
      spin['to'] = 0
      spin['from'] = -100
      spin['increment'] = 0.2
      slide['to'] = 0
      slide['from'] = -100
      slide['digits'] = 3
      slide['resolution'] = 0.2      
      spin.grid(row=3, column=1, sticky='news')
      slide.grid(row=3, column=2, sticky='news')

      # ts
      def update_ts(value):
         self.var_ts.set( float(value) )      
         self.show_mbolmodel()   
      l = Tk.Label(menu,text='T turn (day) =',fg="black")      
      l.grid(row=4,column=0)
      self.var_ts = Tk.DoubleVar(value=60) 
      spin = Tk.Spinbox(menu, textvariable=self.var_ts, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_ts, orient='horizontal',
                       length=200, command=update_ts)
      spin['to'] = 120
      spin['from'] = 20
      spin['increment'] = 1
      slide['to'] = 120
      slide['from'] = 20
      slide['digits'] = 2
      slide['resolution'] = 1
      spin.grid(row=4, column=1, sticky='news')
      slide.grid(row=4, column=2, sticky='news')
      
      # taum
      def update_taum(value):
         self.var_taum.set( float(value) )
         # fix vexp
         mej, ek = radioactivemodels.taum_to_Mej_Ek( self.var_taum.get(), self.var_vexp.get() )      
         self.var_mej.set( mej )
         self.var_ek.set( ek )
         self.var_t0.set( radioactivemodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.show_mbolmodel()
      l = Tk.Label(menu,text='Arnett taum (day) =',fg="black")      
      l.grid(row=5,column=0)
      self.var_taum = Tk.DoubleVar(value=radioactivemodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()))
      spin = Tk.Spinbox(menu, textvariable=self.var_taum, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_taum, orient='horizontal',
                       length=200, command=update_taum)
      spin['to'] = 100
      spin['from'] = .1
      spin['increment'] = 0.1
      slide['to'] = 100
      slide['from'] = .1
      slide['digits'] = 4
      slide['resolution'] = 0.1
      spin.grid(row=5, column=1, sticky='news')
      slide.grid(row=5, column=2, sticky='news')
      
      # t0
      def update_t0(value):
         self.var_t0.set( float(value) )
         # fix vexp
         mej, ek = radioactivemodels.t0_to_Mej_Ek( self.var_t0.get(), self.var_vexp.get() )
         self.var_mej.set( mej )
         self.var_ek.set( ek )
         self.var_taum.set( radioactivemodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()) )
         self.show_mbolmodel()
      l = Tk.Label(menu,text='Tail t0 (day) =',fg="black")      
      l.grid(row=6,column=0)
      self.var_t0 = Tk.DoubleVar(value=radioactivemodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()))
      spin = Tk.Spinbox(menu, textvariable=self.var_t0, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_t0, orient='horizontal',
                       length=200, command=update_t0)
      spin['to'] = 2000
      spin['from'] = 0
      spin['increment'] = 0.1
      slide['to'] = 2000
      slide['from'] = 0
      slide['digits'] = 4
      slide['resolution'] = 0.1
      spin.grid(row=6, column=1, sticky='news')
      slide.grid(row=6, column=2, sticky='news')
      
      # vexp
      def update_vexp(value):
         self.var_vexp.set( float(value) )
         # fix taum
         mej, ek = radioactivemodels.taum_to_Mej_Ek( self.var_taum.get(), self.var_vexp.get() )
         self.var_mej.set( mej )
         self.var_ek.set( ek )
         self.var_t0.set( radioactivemodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.show_mbolmodel()
      l = Tk.Label(menu,text='V phot (10**3 km/s) =',fg="black")      
      l.grid(row=7,column=0)
      self.var_vexp = Tk.DoubleVar(value=radioactivemodels.Mej_Ek_to_vej(self.var_mej.get(), self.var_ek.get()))
      spin = Tk.Spinbox(menu, textvariable=self.var_vexp, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_vexp, orient='horizontal',
                       length=200, command=update_vexp)
      spin['to'] = 30
      spin['from'] = 2
      spin['increment'] = 0.2
      slide['to'] = 30
      slide['from'] = 2
      slide['digits'] = 4
      slide['resolution'] = 0.2
      spin.grid(row=7, column=1, sticky='news')
      slide.grid(row=7, column=2, sticky='news')
      
      # SBO me
      def update_me(value):
         self.var_me.set( float(value) )      
         self.show_mbolmodel()
      l = Tk.Label(menu,text='M E (Msun) =',fg="brown")      
      l.grid(row=8,column=0)
      self.var_me = Tk.DoubleVar(value=.1)
      spin = Tk.Spinbox(menu, textvariable=self.var_me, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_me, orient='horizontal',
                       length=200, command=update_me)
      spin['to'] = 1
      spin['from'] = 0
      spin['increment'] = 0.01
      slide['to'] = 1
      slide['from'] = 0
      slide['digits'] = 4
      slide['resolution'] = 0.01
      spin.grid(row=8, column=1, sticky='news')
      slide.grid(row=8, column=2, sticky='news')
      
      # SBO re
      def update_re(value):
         self.var_re.set( float(value) )      
         self.show_mbolmodel()
      l = Tk.Label(menu,text='R E (Rsun) =',fg="brown")      
      l.grid(row=9,column=0)
      self.var_re = Tk.DoubleVar(value=100)
      spin = Tk.Spinbox(menu, textvariable=self.var_re, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_re, orient='horizontal',
                       length=200, command=update_re)
      spin['to'] = 1000
      spin['from'] = .1
      spin['increment'] = 1
      slide['to'] = 1000
      slide['from'] = .1
      slide['digits'] = 4
      slide['resolution'] = 1
      spin.grid(row=9, column=1, sticky='news')
      slide.grid(row=9, column=2, sticky='news')
      
      # SBO ee
      def update_ee(value):      
         self.var_ee.set( float(value) )      
         self.show_mbolmodel()    
      l = Tk.Label(menu,text='log E E (foe) =',fg="brown")      
      l.grid(row=10,column=0)
      self.var_ee = Tk.DoubleVar(value=np.log10(1e-2))
      spin = Tk.Spinbox(menu, textvariable=self.var_ee, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_ee, orient='horizontal',
                       length=200, command=update_ee)      
      spin['to'] = np.log10(1)
      spin['from'] = np.log10(1e-3)
      spin['increment'] = .1
      slide['to'] = np.log10(1)
      slide['from'] = np.log10(1e-3)
      slide['digits'] = 4
      slide['resolution'] = .1
      spin.grid(row=10, column=1, sticky='news')
      slide.grid(row=10, column=2, sticky='news')
      
      # Arnett/tail
      self.showarnettVar = Tk.IntVar(value=0)            
      showarnettButton = Tk.Checkbutton(menu, text="Arnett model",variable=self.showarnettVar)      
      showarnettButton.grid(row=11,column=0,sticky=Tk.W)
      self.showtailVar = Tk.IntVar(value=0)
      showtailButton = Tk.Checkbutton(menu, text="Tail model",variable=self.showtailVar)
      showtailButton.grid(row=11,column=1,sticky=Tk.W)
      self.showsboVar = Tk.IntVar(value=0)
      showtailButton = Tk.Checkbutton(menu, text="SBO model",variable=self.showsboVar)
      showtailButton.grid(row=11,column=2,sticky=Tk.W)   

   def show_plmodel(self):
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: set object first')
         return
      a = self.var_a1.get()
      texp = self.var_texp1.get()
      alpha = self.var_alpha.get()
      c = self.var_c1.get()     
      showpl = self.showplVar.get()      
      times_pl = np.arange(
         snobject_pars['multiband_early_xrangep'][0],
         snobject_pars['multiband_early_xrangep'][1],
         .1,
      )
      fm, efm = self.ztfp.data[objid]._flux_at('r', 0, interpolation=None)            
      if fm is None: fm = __lc['flux'].max()
      self.Message('A=%.2f texp=%.2f alpha=%.2f C=%.2f'% (a,texp,alpha,c))
      
      # plot
      if not 'wref' in self.__dict__: self.wref = []
      else:
         for w in self.wref:
            try: self.ax.lines.remove(w())
            except: pass
      if showpl == 1:         
         flux = plfunc(times_pl, texp, a, alpha, c)
         
         self.ax.plot(times_pl + self.ztfp.data[objid].t0 - 2458000,
                      flux/fm*snobject_pars['flux_scale'], color='r', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) ) 
         self.canvas.draw_idle()
         
   def show_bazinmodel(self):
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: set object first')
         return
      a = self.var_a.get()
      trise = self.var_trise.get()
      tfall = self.var_tfall.get()
      c = self.var_c.get()
      dt = self.var_dt.get()
      showbazin = self.showbazinVar.get()      
      times_bazin = np.arange(
         snobject_pars['multiband_main_xrangep'][0],
         snobject_pars['multiband_main_xrangep'][1],
         .1,
      )
      fm, efm = self.ztfp.data[objid]._flux_at('r', 0, interpolation=None)            
      if fm is None: fm = __lc['flux'].max()
      self.Message('A=%.2f dt=%.2f trise=%.2f tfall=%.2f C=%.2f'% (a,dt,trise,tfall,c))
      
      # plot
      if not 'wref' in self.__dict__: self.wref = []
      else:
         for w in self.wref:
            try: self.ax.lines.remove(w())
            except: pass
      if showbazin == 1:         
         flux = bazinfunc(times_bazin, a, dt, tfall, trise, c)         
            
         self.ax.plot(times_bazin + self.ztfp.data[objid].t0 - 2458000,
                      flux/fm*snobject_pars['flux_scale'], color='r', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) ) 
         self.canvas.draw_idle()
         
   def show_mbolmodel(self):
      mni = self.var_mni.get()
      mej = self.var_mej.get()
      ek = self.var_ek.get()
      texp = self.var_texp.get()
      ts = self.var_ts.get()
      taum = self.var_taum.get()
      t0 = self.var_t0.get()
      vexp = self.var_vexp.get()
      me = self.var_me.get()
      re = self.var_re.get()
      ee = 10**(self.var_ee.get())
      showarnett = self.showarnettVar.get()
      showtail = self.showtailVar.get()
      showsbo = self.showsboVar.get()
      times_sbo = np.arange(0, 10, 1e-1)
      times_arnett = np.arange(0, ts, 1)
      times_tail = np.arange(ts, 120, 1)
      self.Message('Mni=%.2f Mej=%.2f Ek=%.2f\n\nvexp=%.2f taum=%.1f t0=%.1f\n\ntexp=%.2f ts=%.2f\n\nMe=%.2f Re=%.1f Ee=%.3f'% (mni,mej,ek,vexp,taum,t0,texp,ts,me,re,ee))
      
      # plot
      if not 'wref' in self.__dict__: self.wref = []
      else:
         for w in self.wref:
            try: self.ax.lines.remove(w())
            except: pass
      if showarnett == 1:
         Larnett = radioactivemodels.Arnett_fit_taum(times_arnett, mni, taum)         
         self.ax.plot(times_arnett+texp, Larnett, color='r', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )
      if showtail == 1:
         Ltail = radioactivemodels.tail_fit_t0(times_tail, mni, t0)
         self.ax.plot(times_tail+texp, Ltail, color='orange', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )
      if showsbo == 1:
         Lsbo = sbomodels.shock_fit(times_sbo, me, re, ee)
         self.ax.plot(times_sbo+texp, Lsbo, color='blue', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )
      if showarnett+showtail+showsbo>0:
         self.canvas.draw_idle()
   
   def Message(self,message,disable=True,restart=True):
      
      self.textbox.config(state='normal',foreground='black',background='SystemButtonFace')
      if restart: self.textbox.delete(1.0,Tk.END)
      self.textbox.insert(Tk.END,message)
      self.msgwindow.config(command=self.textbox.yview)
      if disable: self.textbox.config(state='disabled')
      
   def quit(self):
      if messagebox.askokcancel("Quit", "Do you want to quit?"):
         self.parent.quit()
         self.parent.destroy()
      
################################   call to main #########################
def main():
   root = Tk.Tk()
   screen_width = root.winfo_screenwidth()
   screen_height = root.winfo_screenheight()
   root.geometry("%dx%d+10+10"%(screen_width/1.5,screen_height/1.2))
   root.configure(background='SystemButtonFace')
   root.title('HAFFET GUI (version '+__version__+')'+' '*5+'@SY')   
   # cannot be resized
   root.resizable(0, 0)
   # layout on the root window
   root.columnconfigure(0, weight=1)
   root.columnconfigure(1, weight=4)
   root.rowconfigure(0, weight=4)
   root.rowconfigure(1, weight=1)  
   default_font = font.nametofont("TkFixedFont")
   default_font.configure(weight="bold")
   root.option_add( "*font", default_font)
   root.protocol('WM_DELETE_WINDOW', root.quit)
   
   app = MainMenu(root, screen_width, screen_height)
   app.mainloop()
