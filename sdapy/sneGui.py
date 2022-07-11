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
from tabulate import tabulate
from weakref import ref
from astropy.time import Time

from sdapy import __version__, snerun, \
   image_tool, constants, read_default
from sdapy.pbar import get_progress_bar
from sdapy.model_fitters import get_pars
from sdapy.models.arnett import functions as arnettmodels
from sdapy.models.tail import functions as tailmodels
from sdapy.models.sbo import functions as sbomodels
from sdapy.models.bazin.functions import bazin as bazinfunc

import matplotlib.pyplot as plt
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
if len(snerun.snelist_pars['idkey1']) > 0:
   searchnames = [snerun.snelist_pars['idkey'], snerun.snelist_pars['idkey1']]
else:
   searchnames = [snerun.snelist_pars['idkey']]

dataoptions = {   
   'Data preparation': 'down/re-load ZTF/ATLAS/externel forced/alert photometry/spectra',
   'Multiband LC fits': 'Fit multiband LCs, e.g. via Bazin fits for main peak, power law for early LCs, etc, with Gaussian Process (george), Monte Carlo (emcee), or minimize (scipy)',
   'Colours/Bolometric': 'calculate magnitude differences bwtween different bands/phases.\n\nestimate host ebv by comparing to intrinstic colors\n\nset peak/explosion epochs.\n\ncalculate bolometric LCs with Blackbody fits, or Lyman 2016 analytic functions',
   'Bolometric LC fits': 'Fit bolometric LCs, e.g. shock cooling breakout fits for early tail, Arnett fits for the main peak, tail fits, etc.',
   'Spectral Line fits': 'Fit for spectral element line emission/absorption, e.g. via Gaussian or viglot',
}

plotoptions = {
   'Flux LC': 'Multi-band Lightcurves in uJy',
   'Mag LC' : 'Multi-band Lightcurves in AB magnitude',   
   'colour': 'Colour evlotion plots',
   'luminosity' : 'Bolometric Lightcurves in erg/s or magnitude',
   'spectra' : 'Spectral full range or specific line range',   
   'separator1': None,
   'SNe scatter': 'Population for 1 parameter (Histogram) or 2 parameters (scatter), show their correlation',
   'separator2': None,
   'Finder': 'create a finder with online sources, e.g. PS1',
   'LC baseline': 'Baseline check plots for ZTF forced photometry',
}

distparameters = [
   'z', 'ra', 'dec', 'mkwebv', 'hostebv', 'sntype', 'dm', 't0', 'texp', 'mpeak',
   'deltam15', 'deltam-10', 'trise', 'Lpeak', 'Mni', 'Ekin', 'Mej', 'vexp',
   'taum', 't0',
]

allparameters = [
   'lc', 'colour', 'mbol'
]

# ztfquery source path
LOCALSOURCE = os.getenv('ZTFDATA',"./Data/")

###############################     MAIN WINDOW ##########################
class MainMenu(Tk.Frame):

   def __init__(self, parent, width, height):      
      
      Tk.Frame.__init__(self, parent)
      self.parent = parent
      self.width = width
      self.height = height
      
      # init sdapy class      
      self.obj = Tk.StringVar()
      self.obj.set(None)      
      self.ztfp = snerun.snelist(print_func=self.Message)
      
      # initialize
      self.ztfp.parse_meta()
      self.ztfp.parse_params()
      
      # init GUI
      self.initUI()
      
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
      tb.grid(row=0,column=0,rowspan=3)
      # meta button
      checkMetaBut = Tk.Button(labelFrame,text="Meta",command=self.showmeta)
      checkMetaBut.grid(row=0, column=1, sticky=Tk.NW)
      # load cache button
      lcBut = Tk.Button(labelFrame,text="load",\
               fg="green", command=self.loadall, activebackground='slate grey')
      lcBut.grid(row=1, column=1, sticky=Tk.NW)
      # save cache button
      scBut = Tk.Button(labelFrame,text="save",\
               fg="orange", command=self.saveall, activebackground='slate grey')
      scBut.grid(row=2, column=1, sticky=Tk.NW)
      
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
      searchFrameBut.grid(row=0, column=1, sticky=Tk.SE)
      
      # search txt
      self.searchFrameTxt = Tk.Entry(objlistFrame, justify=Tk.LEFT, width=12)
      self.searchFrameTxt.grid(row=1, column=0, sticky=Tk.SW)      
      
      # add/del
      inputFrameBut = Tk.Button(objlistFrame,text="Add",command=self.addobj)
      inputFrameBut.grid(row=1, column=1, sticky=Tk.SW)
      inputFrameBut = Tk.Button(objlistFrame,text="Del",command=self.delobj)
      inputFrameBut.grid(row=2, column=1, sticky=Tk.SW)      
      
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

      self.fitVar = Tk.IntVar(value=0)
      fitButton = Tk.Checkbutton(configFrame, text="Fit",variable=self.fitVar)      
      fitButton.grid(row=0,column=1,sticky=Tk.W)
      
      self.markerVar = Tk.IntVar(value=0)
      markerButton = Tk.Checkbutton(configFrame, text="marker",variable=self.markerVar)
      markerButton.grid(row=1,column=0,sticky=Tk.W)

      self.zoomVar = Tk.IntVar(value=0)
      zoomButton = Tk.Checkbutton(configFrame, text="zoom",variable=self.zoomVar)
      zoomButton.grid(row=1,column=1,sticky=Tk.W)
      
      def showconfig():
         self.Message( 'for general: %s\n\nfor single: %s' % (snerun.snelist_pars, snerun.snobject_pars ))
      defBut = Tk.Button(configFrame,text="default",\
                     fg="blue", command=showconfig, activebackground='slate grey')      
      defBut.grid(row=2, column=0, sticky=Tk.SW)
      
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
         elif hasattr(artist,'get_xy'): x,y = map(lambda v:np.append([],v),artist.get_xy().T)

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
         if self.currentfigure in ['Mag LC', 'contour']: # invert y axis 
            ax.set_ylim(yy,y0)
         else:
            ax.set_ylim(y0,yy)         
      self.canvas.draw()
      
   def motion(self, e):
      x, y = e.xdata, e.ydata      
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
               if self.currentfigure in ['Mag LC', 'contour']: # invert y axis
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
   def showmeta(self, objid=None):
      if objid is None:
         self.Message( tabulate(self.ztfp.meta, headers = 'keys', tablefmt = 'grid') )
      else:
         self.Message( tabulate(self.ztfp.meta.query('%s==@objid'%(snerun.snelist_pars['idkey'])),
                                headers = 'keys', tablefmt = 'grid') )      
      
   def loadall(self):
      objlist = []      
      with get_progress_bar(True, len(self.ztfp.meta.index)) as pbar:            
         for i, objid in enumerate(self.ztfp.meta.index):                        
            # load data
            loaded = self.ztfp.load_data(objid)
            if loaded: objlist.append(objid)
            
            # update
            pbar.update(1)
            
      if len(objlist) == 0:
         self.Message('No cached objects found')
         return
      self.Message('%i objects have been reloaded'%len(objlist))
      
      # set objid
      self.setobjid(objlist[0])
      
      # update option menu
      menu = self.objlistMenu["menu"]
      menu.delete(0, "end")
      for t in objlist:
         menu.add_command(label=t,command=lambda k=t: self.setobjid(k))
         
      # renew header infos
      self.metainfo.set('Obj: %s \n(%i of %i)' % (self.obj.get(), len(self.ztfp.data), len(self.ztfp.meta)))
            
   def saveall(self):
      
      objlist = list(self.ztfp.data.keys())

      if self.clobberVar.get()==0: clobber = False
      else: clobber = True
      
      for objid in objlist:
         saved = self.ztfp.save_data(objid, force=clobber)               
   
   def addobj(self):      
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: set object first')
         return
      objlist = list(self.ztfp.data.keys())
      if objid in objlist and self.clobberVar.get()==0:
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
      for t in objlist:
         menu.add_command(label=t,command=lambda k=t: self.setobjid(k))
         
      # renew header infos
      self.metainfo.set('Obj: %s \n(%i of %i)' % (self.obj.get(), len(self.ztfp.data), len(self.ztfp.meta)))
      
   def delobj(self):
      objid = self.obj.get()
      if objid in ['','None']:
         self.Message('Warning: No objects left')
         return
      del self.ztfp.data[objid]      
      objlist = list(self.ztfp.data.keys())
      self.Message('deleted %s, %i objects left '%(objid, len(objlist)))
      if len(objlist)>0: self.obj.set(objlist[0])
      else: self.obj.set(None)
      # update option menu
      menu = self.objlistMenu["menu"]
      menu.delete(0, "end")
      for t in objlist:
         menu.add_command(label=t,command=lambda k=t: self.setobjid(k))      
      # renew header infos
      self.metainfo.set('Obj: %s \n(%i of %i)' % (self.obj.get(), len(self.ztfp.data), len(self.ztfp.meta)))
         
   def setobjid(self, obj):
      self.obj.set(obj)
      # renew header info
      self.metainfo.set('Obj: %s \n(%i of %i)' % (self.obj.get(), len(self.ztfp.data), len(self.ztfp.meta)))
      
   def SearchObj(self):
      _str = self.searchFrameTxt.get()
      if _str == '':  # empty for tutorial
         self.Message('input %s' % searchnames)
         return
      objid = self.matchstr(_str, True)
      if objid is not None:
         self.obj.set(objid)
         
   def matchstr(self, _str, verbose=False):
      if len(self.ztfp.meta) == 0:
         msg = 'Warning:\n\nno meta, parse it or manually input parameters'
         flag = None
      else:
         candlist = []
         for i, objid in enumerate(self.ztfp.meta.index):
            for seachname in searchnames:
               if not seachname in self.ztfp.meta.keys(): continue
               kid = self.ztfp.meta[seachname][i]
               if _str in objid or _str in kid: candlist.append(objid)               
         if len(candlist) == 0:
            msg = 'Warning:\n\nno objects matched to %s'%_str
            flag = None
         elif len(candlist) == 1:
            msg = 'Nice:\n\n%s was matched successfully as %s\n\nInfo:\n\n%s'%\
               (_str,candlist[0], tabulate(self.ztfp.meta.query('%s=="%s"'%(snerun.snelist_pars['idkey'], candlist[0])), headers = 'keys', tablefmt = 'grid'))
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
      if self.optvar.get() == 'Colours/Bolometric':
         if not 'bbtop' in self.__dict__: self.bb_popup(objid)
         elif self.bbtop.winfo_exists() == 0: self.bb_popup(objid)
         else: self.Message('Parse Colours/Bolometric topup window already exists')
      if self.optvar.get() == 'Bolometric LC fits':
         if not 'bbtop1' in self.__dict__: self.bb_popup1(objid)
         elif self.bbtop1.winfo_exists() == 0: self.bb_popup1(objid)
         else: self.Message('Fit Bolometric topup window already exists')          
      if self.optvar.get() == 'Spectral Line fits':
         if not 'spectop' in self.__dict__: self.spectrafit_popup(objid)
         elif self.spectop.winfo_exists() == 0: self.spectrafit_popup(objid)
         else: self.Message('Spectral Line Fit topup window already exists')           
         
      self.Message('task %s complete' % self.optvar.get())
      # reset action var
      self.optvar.set(None)
      
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
      if self.plotvar.get() == 'Flux LC':
         if not 'lc' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse lc for %s first'%objid)
            return
         ax = self.ztfp.data[objid].ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax(
            show_title=False, show_legend=True, ylabel_2right=False,
            x0=2458000, source=None
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show Flux lcs for %s'%objid)

         if self.fitVar.get() == 1:
            if not 'fluxtop' in self.__dict__: self.fluxfit_popup()
            elif self.fluxtop.winfo_exists() == 0: self.fluxfit_popup()
            else: self.Message('Flux plot window already exists')
         
      if self.plotvar.get() == 'Mag LC':
         if not 'lc' in self.ztfp.data[objid].__dict__.keys():
            self.Message('Warning: parse lc for %s first'%objid)
            return          
         ax = self.ztfp.data[objid].ax2 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax2(
            show_title=False, show_legend=True, ylabel_2right=False,
            corr_mkw=False, corr_host=False, source=None
            )                                    
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show Mag lcs for %s'%objid)

      if self.plotvar.get() == 'colour': 
         ax = self.ztfp.data[objid].ax3 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax3(show_title=False, show_legend=True, ylabel_2right=False,
                                    corr_mkw=True, corr_host=False, source=None) 
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
         self.ztfp.data[objid]._ax4(show_title=False, show_legend=True,
                                    ylabel_2right=False, logscale=True,)
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
         baseline = self.ztfp.data[objid].calibrate_baseline(ax=ax, key='fcqfid', source='ztffp',
                                       xmin=-100, xmax=-20, ax_xlim=None, ax_ylim=None)
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('calibrate baseline for %s %s' % (objid,baseline))

      
      if self.plotvar.get() == 'spectra':                   
         ax = self.ztfp.data[objid].ax1 = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax1(
            show_title=False, show_legend=True, source=None,
         )
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
      self.plotvar.set(None)
      
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
      if 'bolfit' in checkitems:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            msg += 'Fits: \n'
            for engine in self.ztfp.data[objid].fitcls:
               if engine not in ['bol_main', 'bol_early', 'bol_late']: continue
               msg += ' - %s\n' % (engine)
               for source in self.ztfp.data[objid].fitcls[engine]:
                  if source == 'mbol':
                     msg += '    Lyman Bol LC \n'
                  else:
                     msg += '    BB Bol LC \n' % model 
                  for model in self.ztfp.data[objid].fitcls[engine][source]:                  
                     msg += '        %s \n' % model
         else:
            msg += 'Fits: - \n'
      if 'sed' in checkitems:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            msg += 'Fits: \n'
            for engine in self.ztfp.data[objid].fitcls:
               if engine not in ['specline']: continue
               msg += ' - %s\n' % (engine)
               for source in self.ztfp.data[objid].fitcls[engine]:                  
                  msg += '    %s \n'%source
         else:
            msg += 'Fits: - \n'
      return msg

   ### popup for scatter ###
   def dist_popup(self):
      if self.clobberVar.get()==0: clobber = False
      else: clobber = True
      
      # pop up window
      self.disttop = Tk.Toplevel(self.parent)
      self.disttop.geometry("500x500+50+50")      
      self.disttop.title('Distribute features for %i objects'%len(self.ztfp.data))

      # 1d
      Frame = Tk.LabelFrame(self.disttop, text="1D histogram")
      Frame.grid(row=0,column=0,sticky=Tk.W)
                  
      _Frame = Tk.LabelFrame(Frame, text="Parameter")
      _Frame.grid(row=0,column=0)
      self.pval = Tk.StringVar()
      self.pval.set( distparameters[0] )
      Menu = Tk.OptionMenu(_Frame, self.pval, *distparameters)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # 2d
      Frame = Tk.LabelFrame(self.disttop, text="2D scatter")
      Frame.grid(row=1,column=0,sticky=Tk.W)
                  
      _Frame = Tk.LabelFrame(Frame, text="Parameter x")
      _Frame.grid(row=0,column=0)
      self.pval1 = Tk.StringVar()
      self.pval1.set( distparameters[0] )
      Menu = Tk.OptionMenu(_Frame, self.pval1, *distparameters)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      _Frame = Tk.LabelFrame(Frame, text="Parameter y")
      _Frame.grid(row=1,column=0)
      self.pval2 = Tk.StringVar()
      self.pval2.set( distparameters[0] )
      Menu = Tk.OptionMenu(_Frame, self.pval2, *distparameters)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)  

      # plot all
      Frame = Tk.LabelFrame(self.disttop, text="Show all")
      Frame.grid(row=2,column=0,sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Type")
      _Frame.grid(row=0,column=0)
      self.alltypeval = Tk.StringVar()
      self.alltypeval.set( allparameters[0] )
      Menu = Tk.OptionMenu(_Frame, self.alltypeval, *allparameters)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # plots
      def run1d():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         p = self.pval.get()        
         self.currentfigure = 'scatter'
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')        
         self.ztfp.show1d(p, syntax=None, nbin=10, fontsize=12, labelpad=12)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, p)
         self.set_limits( ax )
         self.Message('1d distribution plots')         
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=run1d)
      inputFrameBut.grid(row=0, column=1, sticky=Tk.W)

      def run2d():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         p1 = self.pval1.get()
         p2 = self.pval2.get()         
         self.currentfigure = 'scatter'
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')        
         self.ztfp.show2d(p1, p2, syntax=None, fontsize=12, labelpad=12)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, '%s vs %s'%(p1,p2))
         self.set_limits( ax )
         self.Message('2d distribution plots')
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=run2d)
      inputFrameBut.grid(row=1, column=1, sticky=Tk.W)

      def runall():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         p = self.alltypeval.get()
         self.currentfigure = 'allplots'
         ax = self.ztfp.ax = self.add_subplot(111, zoomable='horizontal', cursor='both')        
         self.ztfp.showall(p, syntax=None, fontsize=12, labelpad=12)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, 'show all %s'%p)
         self.set_limits( ax )
         self.Message('show %s plots for all objects'%p)
      inputFrameBut = Tk.Button(self.disttop, text="Run", command=runall)
      inputFrameBut.grid(row=2, column=1, sticky=Tk.W)
      
   ### popup for Data ###
   def data_popup(self, objid):
      if self.clobberVar.get()==0: clobber = False
      else: clobber = True
      
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
      
      def showauth(): self.Message( snerun.keypairs )
      def showmeta():
         msg = ''
         for k in ['objid', 'z', 'ra', 'dec', 'mkwebv', 'hostebv',
                   'sntype', 'dm', 't0', 'texp', 'fpeak']:
            msg += '%s: %s \n' % (k, str(self.ztfp.data[objid].__dict__[k]))
         self.Message( msg )
      authbutton = Tk.Button(row11,text="auth",command=showauth)
      authbutton.grid(row=0, column=0)
      metabutton = Tk.Button(row11,text="Info",command=lambda k=objid: self.showmeta(k))
      metabutton.grid(row=0, column=1)
      metabutton = Tk.Button(row11,text="meta",command=showmeta)
      metabutton.grid(row=0, column=2)
      
      if len(snerun.snobject_pars['mjdstart']) > 0:
         mjdstart = '%.1f'%(snerun.snobject_pars['mjdstart'])
      elif self.ztfp.data[objid].t0 > 2400000:
         mjdstart = '%.1f' %(float(self.ztfp.data[objid].t0) - 2400000.5 + float(snerun.snobject_pars['dstart']))
      else:
         mjdstart = None
      if len(snerun.snobject_pars['mjdend']) > 0:
         mjdend = '%.1f'%snerun.snobject_pars['mjdend']
      elif self.ztfp.data[objid].t0 > 2400000:
         mjdend = '%.1f'%(float(self.ztfp.data[objid].t0) - 2400000.5 + float(snerun.snobject_pars['dend']))
      else:
         mjdend = None
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
      
      snrt = snerun.snobject_pars['snrt']
      l = Tk.Label(row11,text='SNR:',fg="black")      
      l.grid(row=2,column=0)
      self.snrentry = Tk.Entry(row11, justify=Tk.LEFT, width=5)
      self.snrentry.insert(0, snrt)
      self.snrentry.grid(row=2, column=1, sticky=Tk.SW)
      l = Tk.Label(row11,text='sigma',fg="black")      
      l.grid(row=2,column=2)
      
      tdbin = snerun.snobject_pars['tdbin']
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
         self.ztfp.data[objid].get_fp_atlas(binDays=int(self.tdbentry.get()), clobber=clobber)
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")
      def downloadatlasforce():
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
         self.ztfp.data[objid].get_alert_ztf(source='marshal')
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
         self.ztfp.data[objid].get_alert_ztf(source='fritz')
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
      # GM spectra
      def loadgmspec():
         self.ztfp.data[objid].get_local_spectra(source='marshal')
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")
      def downloadgmspec():
         self.ztfp.data[objid].query_spectra(source='marshal')
      data1 = Tk.Label(self.datatop,text='ZTF Growth marshal Spectra',fg="black")
      data1.grid(row=5,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadgmspec)
      data1button.grid(row=5, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadgmspec)
      data2button.grid(row=5, column=2)
      # fritz spectra
      def loadfritzspec():
         self.ztfp.data[objid].get_local_spectra(source='fritz')
         self.datainfo1.set(self.checkdata(objid, ['data']))
         tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")
      def downloadfritzspec():
         self.ztfp.data[objid].query_spectra(source='fritz')
      data1 = Tk.Label(self.datatop,text='ZTF fritz Spectra',fg="black")
      data1.grid(row=6,column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadfritzspec)
      data1button.grid(row=6, column=1)
      data2button = Tk.Button(self.datatop,text="download",command=downloadfritzspec)
      data2button.grid(row=6, column=2)
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
            self.ztfp.data[objid].get_external_phot(filename, source)
            self.datainfo1.set(self.checkdata(objid, ['data']))
            tb = Tk.Label(self.datatop,textvariable=self.datainfo1,fg="red")             
      data1 = Tk.Label(self.datatop,text='external Photometry',fg="black")
      data1.grid(row=7,column=0)
      dataFrame = Tk.LabelFrame(self.datatop, text="source name")
      dataFrame.grid(row=7,column=1)
      self.extphotentry = Tk.Entry(dataFrame, width=10)      
      self.extphotentry.grid(row=0, column=0)   
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadextphot)
      data1button.grid(row=7, column=2)         
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
      data1.grid(row=8,column=0)
      dataFrame = Tk.LabelFrame(self.datatop, text="epochs")
      dataFrame.grid(row=8,column=1)
      self.extspecentry = Tk.Entry(dataFrame, width=10)
      self.extspecentry.grid(row=0, column=0)
      data1button = Tk.Button(self.datatop,text="loadlocal",command=loadextspec)
      data1button.grid(row=8, column=2)
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
      data1.grid(row=9,column=0)
      self.inputFrameTxt = Tk.Entry(self.datatop, justify=Tk.LEFT, width=12)
      self.inputFrameTxt.grid(row=9, column=1, sticky=Tk.SW)
      inputFrameBut = Tk.Button(self.datatop,text="submit",command=InputMeta)
      inputFrameBut.grid(row=9, column=2, sticky=Tk.SE)

      def cliplc():  self.ztfp.data[objid].clip_lc()
      inputFrameBut = Tk.Button(self.datatop,text="clip LC",command=cliplc)
      inputFrameBut.grid(row=10, column=0, sticky=Tk.SW)

      def rapid():  self.ztfp.data[objid].rapid()
      inputFrameBut = Tk.Button(self.datatop,text="Guess SNtype",command=rapid)
      inputFrameBut.grid(row=10, column=0, sticky=Tk.SW)
      
      def midtexp(): self.ztfp.data[objid].set_texp_midway()
      inputFrameBut = Tk.Button(self.datatop,text="set texp with (non) detecitons",command=midtexp)
      inputFrameBut.grid(row=11, column=0, sticky=Tk.SW)
      
   ### popup for Data 1 ###
   def data_popup1(self, objid):      
      # pop up window
      self.datatop1 = Tk.Toplevel(self.parent)
      self.datatop1.geometry("500x500+50+50")      
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
      
      nsteps = snerun.snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=1,column=0)
      self.nstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepsentry.insert(0, nsteps)
      self.nstepsentry.grid(row=1, column=1)
      
      nsteps_burnin = snerun.snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=2,column=0)
      self.nstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepburnsentry.insert(0, nsteps_burnin)
      self.nstepburnsentry.grid(row=2, column=1)

      nwalkers = snerun.snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=1,column=2)
      self.nwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nwalkerentry.insert(0, nwalkers)
      self.nwalkerentry.grid(row=1, column=3)
      
      if 'lc' in self.ztfp.data[objid].__dict__:
         colours = list(np.unique(self.ztfp.data[objid].lc['filter']))
      else:
         colours = snerun.snobject_pars['plot_bands']
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
      button = Tk.Button(Frame,text="u",command=updateinfos)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=1, sticky=Tk.SW)
      
      # Gaussian process
      Frame = Tk.LabelFrame(self.datatop1, text="Gaussian Process")
      Frame.grid(row=1,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snerun.snobject_pars['gp_fitr']      
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
      
      gpmin, gpmax = snerun.snobject_pars['gp_plotr']      
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
      self.gproutine.set( snerun.snobject_pars['gp_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.gproutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Kernel")
      _Frame.grid(row=1,column=3)
      _list = ['matern52', 'matern32', 'squaredexp']
      self.kernel = Tk.StringVar()
      self.kernel.set( snerun.snobject_pars['kernel'] )
      Menu = Tk.OptionMenu(_Frame, self.kernel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Mean")
      _Frame.grid(row=0,column=4)
      _list = ['mean', 'gaussian', 'bazin', 'villar']
      self.gpmean = Tk.StringVar()
      self.gpmean.set( snerun.snobject_pars['gp_mean'] )
      Menu = Tk.OptionMenu(_Frame, self.gpmean, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)   
      
      _Frame = Tk.LabelFrame(Frame, text="Fix scale")
      _Frame.grid(row=1,column=4)
      _list = [False, True]
      self.gpfixscale = Tk.StringVar()
      self.gpfixscale.set( snerun.snobject_pars['fix_scale'] )
      Menu = Tk.OptionMenu(_Frame, self.gpfixscale, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)   
      
      # Power law fits
      Frame = Tk.LabelFrame(self.datatop1, text="Multiband early")
      Frame.grid(row=2,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snerun.snobject_pars['multiband_early_xrange']      
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
      
      gpmin, gpmax = snerun.snobject_pars['multiband_early_xrangep']      
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
      self.plroutine.set( snerun.snobject_pars['multiband_early_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.plroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = ['pl', 'pl2']
      self.plmodel = Tk.StringVar()
      self.plmodel.set( 'pl' )
      Menu = Tk.OptionMenu(_Frame, self.plmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # Bazin fits
      Frame = Tk.LabelFrame(self.datatop1, text="Multiband main")
      Frame.grid(row=3,column=0,columnspan=2,sticky=Tk.W)

      gpmin, gpmax = snerun.snobject_pars['multiband_main_xrange']      
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
      
      gpmin, gpmax = snerun.snobject_pars['multiband_main_xrangep']      
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
      self.bzroutine.set( snerun.snobject_pars['multiband_main_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.bzroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = ['bz', 'bz2']
      self.bzmodel = Tk.StringVar()
      self.bzmodel.set( 'bz' )
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
         self.ztfp.data[objid].run_gp(
            gp_bands = gp_bands,
            gp_fitr = gp_fitr,                                      
            kernel = self.kernel.get(),
            fix_scale = self.gpfixscale.get(),
            gp_mean = self.gpmean.get(),
            opt_routine = self.gproutine.get(),
            nwalkers = int(self.nwalkerentry.get()),
            nsteps = int(self.nstepsentry.get()),
            nsteps_burnin = int(self.nstepburnsentry.get()),
            gp_redo = clobber,
            source = source,
         )
         self.datainfo2.set(self.checkdata(objid, ['fits']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo2, fg="red", justify= Tk.LEFT)
      button = Tk.Button(self.datatop1,text="run GP",command=rungp)
      button.grid(row=4, column=0, sticky=Tk.W)
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
      button1 = Tk.Button(self.datatop1,text="run early",command=runpl)
      button1.grid(row=5, column=0, sticky=Tk.W)
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
      button2 = Tk.Button(self.datatop1,text="run main",command=runbz)
      button2.grid(row=6, column=0, sticky=Tk.W)

      # 
      def gpt0(): self.ztfp.data[objid].set_peak_gp(snerun.snobject_pars['set_tpeak_filter'])
      inputFrameBut = Tk.Button(self.datatop1, text="set t0", command=gpt0)
      inputFrameBut.grid(row=4, column=1, sticky=Tk.E)

      def maint0():
         self.ztfp.data[objid].set_peak_multiband_main(snerun.snobject_pars['set_tpeak_filter'])
      inputFrameBut = Tk.Button(self.datatop1, text="set t0", command=maint0)
      inputFrameBut.grid(row=6, column=1, sticky=Tk.E)      

      def pltexp():
         self.ztfp.data[objid].set_texp_pl(snerun.snobject_pars['set_texp_filter'])
      inputFrameBut = Tk.Button(self.datatop1, text="set texp", command=pltexp)
      inputFrameBut.grid(row=5, column=1, sticky=Tk.E)

      # corner plots
      def gpcontours():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, gp=True)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)         
      inputFrameBut = Tk.Button(self.datatop1, text="Contours", command=gpcontours)
      inputFrameBut.grid(row=4, column=1, sticky=Tk.W)

      def maincontours():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, engine='multiband_main', model=None, source=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.datatop1, text="Contours", command=maincontours)
      inputFrameBut.grid(row=6, column=1, sticky=Tk.W)      

      def plcontours():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, engine='multiband_early', model=None, source=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.datatop1, text="Contours", command=plcontours)
      inputFrameBut.grid(row=5, column=1, sticky=Tk.W)
      
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
                  
      row = Tk.LabelFrame(row11, text="Match colours")
      row.grid(row=0,column=0)
      
      self.interplist = Tk.Listbox(row, selectmode = "multiple")  
      self.interplist.grid(row=0,column=0)              
      for each_item in ["bin","gp","fit"]: self.interplist.insert(Tk.END, each_item)
      
      # row1 right
      row11 = Tk.Frame(self.bbtop)
      row11.grid(row=0,column=2,sticky=Tk.W)
      
      tdbin = snerun.snobject_pars['tdbin']
      l = Tk.Label(row11,text='bin:',fg="black")      
      l.grid(row=0,column=0)
      self.tdbentry1 = Tk.Entry(row11, justify=Tk.LEFT, width=5)
      self.tdbentry1.insert(0, tdbin)
      self.tdbentry1.grid(row=0, column=1, sticky=Tk.SW)
      l = Tk.Label(row11,text='days',fg="black")      
      l.grid(row=0,column=2)

      self.corrmkwebv = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(row11, text="corr mkw",variable=self.corrmkwebv) 
      cButton.grid(row=1,column=0,columnspan=3,sticky=Tk.W)
      
      self.corrhostebv = Tk.IntVar(value=0)
      cButton = Tk.Checkbutton(row11, text="corr host",variable=self.corrhostebv)      
      cButton.grid(row=2,column=0,columnspan=3,sticky=Tk.W)
      
      # calculate colours
      if snerun.snobject_pars['color_bands'] is None:
         colours = ['g','r']
      else:
         colours = snerun.snobject_pars['color_bands']              
      data1 = Tk.Label(self.bbtop,text='calculate colors',fg="black")
      data1.grid(row=1,column=0)        
      colourFrame = Tk.LabelFrame(self.bbtop, text="2 bands")
      colourFrame.grid(row=1,column=1)
      self.colorentry1 = Tk.Entry(colourFrame, width=10)
      self.colorentry1.insert(0, colours)
      self.colorentry1.grid(row=0, column=0)      
      def makecolor():
         bands = self.colorentry1.get().split()
         cinterp = []
         for index in self.interplist.curselection():
            if self.interplist.get(index) == 'bin': cinterp.append(1)
            if self.interplist.get(index) == 'fit': cinterp.append(2)
            if self.interplist.get(index) == 'gp':  cinterp.append(3) 
         if len(bands) != 2:
            self.Message('Need/Only two bands needed')
            return
         if len(cinterp) == 0:
            self.Message('Select at least one option to match colours')
            return
         self.ztfp.data[objid].calc_colors(
            color_bands = bands,
            tdbin = int(self.tdbentry1.get()),
            color_interp = cinterp,
         )
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(self.bbtop,textvariable=self.datainfo3,fg="red")
      data1button = Tk.Button(self.bbtop,text="submit",command=makecolor)
      data1button.grid(row=1, column=2)
      
      # make lyman bolometric
      if snerun.snobject_pars['color_bands'] is None:
         colours = ['g','r']
      else:
         colours = snerun.snobject_pars['color_bands'] 
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
         self.ztfp.data[objid].lyman_bol(
            color_bands=bands
         )
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(self.bbtop,textvariable=self.datainfo3,fg="red")         
      data1button = Tk.Button(self.bbtop,text="submit",command=makebc)
      data1button.grid(row=2, column=2)
      
      # make BB bolometric
      if snerun.snobject_pars['bb_bands'] is None:
         colours = ['g','r','i']
      else:
         colours = snerun.snobject_pars['bb_bands']   
      data1 = Tk.Label(self.bbtop,text='Bolometric with BB',fg="black")
      data1.grid(row=3,column=0)        
      bbFrame = Tk.LabelFrame(self.bbtop, text="multi bands")
      bbFrame.grid(row=3,column=1)      
      self.bbentry = Tk.Entry(bbFrame, width=10)
      self.bbentry.insert(0, colours)
      self.bbentry.grid(row=0, column=0)
      
      def makebc1():
         bands = self.bbentry.get().split()
         cinterp = []
         for index in self.interplist.curselection():
            if self.interplist.get(index) == 'bin': cinterp.append(1)
            if self.interplist.get(index) == 'fit': cinterp.append(2)
            if self.interplist.get(index) == 'gp':  cinterp.append(3) 
         if len(bands) <= 2:
            self.Message('Need at least three bands')
            return
         if len(cinterp) == 0:
            self.Message('Select at least one option to match colours')
            return
         self.ztfp.data[objid].bb_colors(
            bb_bands=bands,
            bb_interp=cinterp,
         )
         self.ztfp.data[objid].bb_bol(
            bb_bands=bands,
         )
         self.datainfo3.set(self.checkdata(objid, ['bc']))
         tb = Tk.Label(self.bbtop,textvariable=self.datainfo3,fg="red")
      data1button = Tk.Button(self.bbtop,text="submit",command=makebc1)
      data1button.grid(row=3, column=2)

      # estimate host
      if snerun.snobject_pars['hostebv_bands'] is None:
         colours = ['g','r']
      else:
         colours = snerun.snobject_pars['hostebv_bands']   
      data1 = Tk.Label(self.bbtop,text='Estimate host ebv',fg="black")
      data1.grid(row=4,column=0)
      hostFrame = Tk.LabelFrame(self.bbtop, text="2 bands")
      hostFrame.grid(row=4,column=1)      
      self.hostentry = Tk.Entry(hostFrame, width=10)
      self.hostentry.insert(0, colours)
      self.hostentry.grid(row=0, column=0)
      
      def esthost():
         bands = self.hostentry.get().split()
         interpolation = self.interplist.get(Tk.ACTIVE)
         if len(bands) != 2:
            self.Message('Need/Only two bands needed')
            return
         hostebv = self.ztfp.data[objid].hostebv
         self.ztfp.data[objid].est_hostebv_with_c10(
            interpolation=interpolation,
            hostebv_bands=bands,
         )
         self.Message( 'Host EBV: %.2f -> %.2f' % (hostebv, self.ztfp.data[objid].hostebv) )
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
         m, em = self.ztfp.data[objid]._mag_at(
            filt, phase,
            interpolation=self.interplist.get(Tk.ACTIVE),
            corr_mkw=corr_mkw,
            corr_host=corr_host,
         )
         self.Message( 'Mag = %s (%s)' % (m, em) )
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
         m, em = self.ztfp.data[objid]._absmag_at(
            filt, phase,
            interpolation=self.interplist.get(Tk.ACTIVE),
            corr_mkw=corr_mkw,
            corr_host=corr_host,
         )
         self.Message( 'Mag = %s (%s)' % (m, em) )
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
         m, em = self.ztfp.data[objid]._rate_at(
            filt, float(phases[0]), float(phases[1]),
            interpolation=self.interplist.get(Tk.ACTIVE)
         )
         self.Message( 'delta Mag (%s-%s) = %s (%s)' % (phases[0], phases[1], m, em) )
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
         m, em = self.ztfp.data[objid]._color_at(
            filters[0], filters[1], float(phase),
            interpolation=self.interplist.get(Tk.ACTIVE),
            corr_mkw=corr_mkw,
            corr_host=corr_host,            
         )
         self.Message( 'delta Mag (%s-%s) = %s (%s)' % (filters[0], filters[1], m, em) )
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
      
      nsteps = snerun.snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=1,column=0)
      self.nstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepsentry.insert(0, nsteps)
      self.nstepsentry.grid(row=1, column=1)
      
      nsteps_burnin = snerun.snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=2,column=0)
      self.nstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepburnsentry.insert(0, nsteps_burnin)
      self.nstepburnsentry.grid(row=2, column=1)

      nwalkers = snerun.snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=1,column=2)
      self.nwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nwalkerentry.insert(0, nwalkers)
      self.nwalkerentry.grid(row=1, column=3)
      
      if snerun.snobject_pars['color_interp'] is None:
         cinterp = [1,2,3]
      else:
         cinterp = snerun.snobject_pars['color_interp']
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
         self.msource.set( sources[0] )         
         # update option menu
         menu = lcsourceMenu["menu"]
         menu.delete(0, "end")
         for t in sources:            
            menu.add_command(label=t,command=lambda k=t: updatesource(k))  
      button = Tk.Button(Frame,text="u",command=updateinfos)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=1, sticky=Tk.SW)
      
      # bol early fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric early (shock cooling)")
      Frame.grid(row=1,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snerun.snobject_pars['bol_early_xrange']      
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
      
      gpmin, gpmax = snerun.snobject_pars['bol_early_xrangep']      
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
      self.sboroutine.set( snerun.snobject_pars['bol_early_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.sboroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = ['shock']
      self.sbomodel = Tk.StringVar()
      self.sbomodel.set( 'shock' )
      Menu = Tk.OptionMenu(_Frame, self.sbomodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      # bol main fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric main (Arnett)")
      Frame.grid(row=2,column=0,columnspan=2,sticky=Tk.W)

      gpmin, gpmax = snerun.snobject_pars['bol_main_xrange']      
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
      
      gpmin, gpmax = snerun.snobject_pars['bol_main_xrangep']      
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
      self.bmroutine.set( snerun.snobject_pars['bol_main_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.bmroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = ['arnett', 'arnett2', 'arnett3', 'arnett4']
      self.bmmodel = Tk.StringVar()
      self.bmmodel.set( 'arnett' )
      Menu = Tk.OptionMenu(_Frame, self.bmmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      # bol late fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric tail (radioactive)")
      Frame.grid(row=3,column=0,columnspan=2,sticky=Tk.W)

      gpmin, gpmax = snerun.snobject_pars['bol_tail_xrange']      
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
      
      gpmin, gpmax = snerun.snobject_pars['bol_tail_xrangep']      
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
      self.tailroutine.set( snerun.snobject_pars['bol_tail_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.tailroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = ['tail', 'tail2', 'tail3', 'tail4']
      self.tailmodel = Tk.StringVar()
      self.tailmodel.set( 'tail' )
      Menu = Tk.OptionMenu(_Frame, self.tailmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      # bol joint fits
      Frame = Tk.LabelFrame(self.bbtop1, text="Bolometric full (Arnett+tail)")
      Frame.grid(row=4,column=0,columnspan=2,sticky=Tk.W)
      
      gpmin, gpmax = snerun.snobject_pars['bol_full_xrange']      
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
      
      gpmin, gpmax = snerun.snobject_pars['bol_full_xrangep']      
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
      self.bfroutine.set( snerun.snobject_pars['bol_full_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.bfroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=1,column=3)
      _list = ['arnett_tail', 'arnett_tail2', 'arnett_tail3', 'arnett_tail4']
      self.bfmodel = Tk.StringVar()
      self.bfmodel.set( 'arnett_tail' )
      Menu = Tk.OptionMenu(_Frame, self.bfmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)

      # done
      def runearlybol():
         return
      button = Tk.Button(self.bbtop1,text="run fits on early phases",command=runearlybol)
      button.grid(row=5, column=0, sticky=Tk.W)
      
      def runmainbol():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.msource.get()
         self.Message('Bol main fits for source: %s'%(source))
         if len(source) == 0 or source == 'None': source=None  
         fit_fitr = [float(self.bmfitminentry.get()),float(self.bmfitmaxentry.get())]
         fit_fitp = [float(self.bmplotminentry.get()),float(self.bmplotmaxentry.get())]
         cinterp = []
         for c in list(self.cinterpentry.get()):
            try: cinterp.append(int(c))
            except: pass
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
      button1 = Tk.Button(self.bbtop1,text="run fits on main peak",command=runmainbol)
      button1.grid(row=6, column=0, sticky=Tk.W)

      def runtailbol():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         fit_fitr = [float(self.tailfitminentry.get()),float(self.tailfitmaxentry.get())]
         fit_fitp = [float(self.tailplotminentry.get()),float(self.tailplotmaxentry.get())]
         cinterp = []
         for c in list(self.cinterpentry.get()):
            try: cinterp.append(int(c))
            except: pass         
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
         )
         self.datainfo4.set(self.checkdata(objid, ['bolfit']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo4,fg="red")      
      button1 = Tk.Button(self.bbtop1,text="run fits on tail",command=runtailbol)
      button1.grid(row=7, column=0, sticky=Tk.W)

      # bol main corner
      def bmaincontours():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, engine='bol_main', model=None, source=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.bbtop1, text="Contours", command=bmaincontours)
      inputFrameBut.grid(row=6, column=1, sticky=Tk.W)

      def btailcontours():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, engine='bol_tail', model=None, source=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)  
      inputFrameBut = Tk.Button(self.bbtop1, text="Contours", command=btailcontours)
      inputFrameBut.grid(row=7, column=1, sticky=Tk.W)
      
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
         sources = np.append(sources, list(self.ztfp.data[objid].spec.data.keys()))
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
      
      #
      nsteps = snerun.snobject_pars['nsteps']
      l = Tk.Label(Frame,text='steps=',fg="black")      
      l.grid(row=1,column=0)
      self.nstepsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepsentry.insert(0, nsteps)
      self.nstepsentry.grid(row=1, column=1)
      
      nsteps_burnin = snerun.snobject_pars['nsteps_burnin']      
      l = Tk.Label(Frame,text='burnin=',fg="black")      
      l.grid(row=2,column=0)
      self.nstepburnsentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nstepburnsentry.insert(0, nsteps_burnin)
      self.nstepburnsentry.grid(row=2, column=1)

      nwalkers = snerun.snobject_pars['nwalkers']      
      l = Tk.Label(Frame,text='walkers=',fg="black")      
      l.grid(row=1,column=2)
      self.nwalkerentry = Tk.Entry(Frame, justify=Tk.LEFT, width=6)
      self.nwalkerentry.insert(0, nwalkers)
      self.nwalkerentry.grid(row=1, column=3)

      # spectra fits
      Frame = Tk.LabelFrame(self.spectop, text="spectral line fits")
      Frame.grid(row=1,column=0,columnspan=2,sticky=Tk.W)

      _Frame = Tk.LabelFrame(Frame, text="line")
      _Frame.grid(row=0,column=0)      
      lines = list(constants.line_location.keys())
      self.specline = Tk.StringVar()
      self.specline.set( lines[0] )
      
      Menu = Tk.OptionMenu(Frame, self.specline, *lines)
      Menu.config(justify=Tk.LEFT,width=4)
      Menu.grid(row=0,column=0,sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Routine")
      _Frame.grid(row=0,column=1)
      _list = ['mcmc', 'minimize', 'leastsq']
      self.specroutine = Tk.StringVar()
      self.specroutine.set( snerun.snobject_pars['bol_early_routine'] )
      Menu = Tk.OptionMenu(_Frame, self.specroutine, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      _Frame = Tk.LabelFrame(Frame, text="Model")
      _Frame.grid(row=0,column=2)
      _list = ['gauss']
      self.specmodel = Tk.StringVar()
      self.specmodel.set( 'gauss' )
      Menu = Tk.OptionMenu(_Frame, self.specmodel, *_list)
      Menu.config(justify=Tk.LEFT,width=8)
      Menu.grid(row=0,column=0, sticky=Tk.W)
      
      # done
      def runfit():
         if self.clobberVar.get()==0: clobber=False
         else: clobber=True
         source = self.specsource.get()
         sn_line = self.specline.get()
         self.Message('%s line fits for source: %s with clobber=%s'%(sn_line,source,clobber))
         if len(source) == 0 or source == 'None': source=None
         if len(sn_line) == 0 or sn_line == 'None': source=None
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
         )
         self.datainfo5.set(self.checkdata(objid, ['sed']))
         tb = Tk.Label(labelFrame,textvariable=self.datainfo5, fg="red", justify= Tk.LEFT)         
      button2 = Tk.Button(self.spectop,text="run fitting",command=runfit)
      button2.grid(row=2, column=0, sticky=Tk.W)
      
      # corner plots
      def speccontours():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         self.currentfigure = 'contour'
         ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid].show_corner(ax, engine='specline', model=None, source=None)
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Contour plots for %s'%objid)
      inputFrameBut = Tk.Button(self.spectop, text="Contours", command=speccontours)
      inputFrameBut.grid(row=2, column=1, sticky=Tk.W)
      
   ### popup for Flux LC ###  
   def fluxfit_popup(self):      
      objid = self.obj.get()
      
      # pop up window
      self.fluxtop = Tk.Toplevel(self.parent)
      self.fluxtop.geometry("500x500+50+50")      
      self.fluxtop.title('Flux model generator')

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
      
      def updatesource(k): self.lcsource1.set( k )
      def updateinfos():
         sources = [None]
         if 'lc' in self.ztfp.data[objid].__dict__:            
            sources = np.append(sources, np.unique(self.ztfp.data[objid].lc['source']))
         self.lcsource1.set( sources[0] )         
         # update option menu
         menu = lcsourceMenu["menu"]
         menu.delete(0, "end")
         for t in sources:            
            menu.add_command(label=t,command=lambda k=t: updatesource(k))  
      button = Tk.Button(self.fluxtop,text="update",command=updateinfos)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=0)

      def updateplot():
         # reset canvas
         try:    self.fig.clear(False)
         except: pass
         source = self.lcsource1.get()
         self.Message('GP for source: %s'%(source))
         if len(source) == 0 or source == 'None': source=None         
         ax = self.ztfp.data[objid].ax = self.add_subplot(111, zoomable='horizontal', cursor='both')
         self.ztfp.data[objid]._ax(
            show_title=False, show_legend=True, ylabel_2right=False,
            x0=2458000, source=source,
         )
         self.currentfiguresize = (ax.get_xlim(), ax.get_ylim())         
         self.set_title(ax, objid)
         self.set_limits( ax )
         self.Message('Show Flux lcs for %s'%objid)
         
      button = Tk.Button(self.fluxtop,text="show",command=updateplot)
      button.config(justify=Tk.LEFT,width=1)
      button.grid(row=0,column=0, sticky=Tk.E)
      
      # Bazin fits
      menu = Tk.LabelFrame(self.fluxtop, text="Bazin model")      
      menu.grid(row=1, column=0) 
      
      # default p
      filt = 'r'
      a_p = constants.a_p0
      dt_p = constants.dt_p0
      trise_p = constants.trise_p0
      tfall_p = constants.tfall_p0
      c_p = constants.c_p0
      if objid is not None:
         if 'fitcls' in self.ztfp.data[objid].__dict__:
            if 'multiband_early' in self.ztfp.data[objid].fitcls:
               for model in self.ztfp.data[objid].fitcls['multiband_early']:
                  try:a_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='a')
                  except:pass
                  try:trise_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='trise')
                  except:pass
                  try:tfall_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='tfall')
                  except:pass
                  try:c_p = self.ztfp.data[objid].fitcls['multiband_early'][model].get_par(filt=filt, parname='c')
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
      spin['to'] = constants.a_bounds[1]
      spin['from'] = constants.a_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.a_bounds[1]
      slide['from'] = constants.a_bounds[0]
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
      spin['to'] = constants.dt_bounds[1]
      spin['from'] = constants.dt_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.dt_bounds[1]
      slide['from'] = constants.dt_bounds[0]
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
      spin['to'] = constants.trise_bounds[1]
      spin['from'] = constants.trise_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.trise_bounds[1]
      slide['from'] = constants.trise_bounds[0]
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
      spin['to'] = constants.tfall_bounds[1]
      spin['from'] = constants.tfall_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.tfall_bounds[1]
      slide['from'] = constants.tfall_bounds[0]
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
      spin['to'] = constants.c_bounds[1]
      spin['from'] = constants.c_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.c_bounds[1]
      slide['from'] = constants.c_bounds[0]
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=4, column=1, sticky='news')
      slide.grid(row=4, column=2, sticky='news')      
      
      # show
      self.showbazinVar = Tk.IntVar(value=0) 
      showbazinButton = Tk.Checkbutton(menu, text="Bazin model",variable=self.showbazinVar)      
      showbazinButton.grid(row=5,column=0,sticky=Tk.W)     
      
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
      self.var_mni = Tk.DoubleVar(value=constants.mni_p0)
      spin = Tk.Spinbox(menu, textvariable=self.var_mni, wrap=True,
                        width=10, command=self.show_mbolmodel)
      def update_mni(value):
         self.var_mni.set( float(value) )      
         self.show_mbolmodel()
      slide = Tk.Scale(menu, variable=self.var_mni, orient='horizontal',
                       length=200, command=update_mni)
      spin['to'] = constants.mni_bounds[1]
      spin['from'] = constants.mni_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.mni_bounds[1]
      slide['from'] = constants.mni_bounds[0]
      slide['digits'] = 3
      slide['resolution'] = 0.01
      spin.grid(row=0, column=1, sticky='news')
      slide.grid(row=0, column=2, sticky='news')

      # mej
      def update_mej(value):
         self.var_mej.set( float(value) )
         self.var_taum.set( arnettmodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()) )
         self.var_t0.set( tailmodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.var_vexp.set( arnettmodels.Mej_Ek_to_vej(self.var_mej.get(), self.var_ek.get()) )     
         self.show_mbolmodel()
      l = Tk.Label(menu,text='M ej (Msun) =',fg="black")      
      l.grid(row=1,column=0)
      self.var_mej = Tk.DoubleVar(value=constants.mej_p0)
      spin = Tk.Spinbox(menu, textvariable=self.var_mej, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_mej, orient='horizontal',
                       length=200, command=update_mej)
      spin['to'] = constants.mej_bounds[1]
      spin['from'] = constants.mej_bounds[0]
      spin['increment'] = 0.1
      slide['to'] = constants.mej_bounds[1]
      slide['from'] = constants.mej_bounds[0]
      slide['digits'] = 3
      slide['resolution'] = 0.1      
      spin.grid(row=1, column=1, sticky='news')
      slide.grid(row=1, column=2, sticky='news')
      
      # ek
      def update_ek(value):
         self.var_ek.set( float(value) )      
         self.var_taum.set( arnettmodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()) )
         self.var_t0.set( tailmodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.var_vexp.set( arnettmodels.Mej_Ek_to_vej(self.var_mej.get(), self.var_ek.get()) )     
         self.show_mbolmodel()
      l = Tk.Label(menu,text='E kin (foe) =',fg="black")      
      l.grid(row=2,column=0)
      self.var_ek = Tk.DoubleVar(value=constants.ek_p0)
      spin = Tk.Spinbox(menu, textvariable=self.var_ek, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_ek, orient='horizontal',
                       length=200, command=update_ek)
      spin['to'] = constants.ek_bounds[1]
      spin['from'] = constants.ek_bounds[0]
      spin['increment'] = 0.1
      slide['to'] = constants.ek_bounds[1]
      slide['from'] = constants.ek_bounds[0]
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
      self.var_texp = Tk.DoubleVar(value=constants.texp_p0)
      spin = Tk.Spinbox(menu, textvariable=self.var_texp, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_texp, orient='horizontal',
                       length=200, command=update_texp)
      spin['to'] = constants.texp_bounds[1]
      spin['from'] = constants.texp_bounds[0]
      spin['increment'] = 0.2
      slide['to'] = constants.texp_bounds[1]
      slide['from'] = constants.texp_bounds[0]
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
      self.var_ts = Tk.DoubleVar(value=constants.ts_p0) 
      spin = Tk.Spinbox(menu, textvariable=self.var_ts, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_ts, orient='horizontal',
                       length=200, command=update_ts)
      spin['to'] = constants.ts_bounds[1]
      spin['from'] = constants.ts_bounds[0]
      spin['increment'] = 1
      slide['to'] = constants.ts_bounds[1]
      slide['from'] = constants.ts_bounds[0]
      slide['digits'] = 2
      slide['resolution'] = 1
      spin.grid(row=4, column=1, sticky='news')
      slide.grid(row=4, column=2, sticky='news')
      
      # taum
      def update_taum(value):
         self.var_taum.set( float(value) )
         # fix vexp
         mej, ek = arnettmodels.taum_to_Mej_Ek( self.var_taum.get(), self.var_vexp.get() )      
         self.var_mej.set( mej )
         self.var_ek.set( ek )
         self.var_t0.set( tailmodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.show_mbolmodel()
      l = Tk.Label(menu,text='Arnett taum (day) =',fg="black")      
      l.grid(row=5,column=0)
      self.var_taum = Tk.DoubleVar(value=arnettmodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()))
      spin = Tk.Spinbox(menu, textvariable=self.var_taum, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_taum, orient='horizontal',
                       length=200, command=update_taum)
      spin['to'] = constants.taum_bounds[1]
      spin['from'] = constants.taum_bounds[0]
      spin['increment'] = 0.1
      slide['to'] = constants.taum_bounds[1]
      slide['from'] = constants.taum_bounds[0]
      slide['digits'] = 4
      slide['resolution'] = 0.1
      spin.grid(row=5, column=1, sticky='news')
      slide.grid(row=5, column=2, sticky='news')
      
      # t0
      def update_t0(value):
         self.var_t0.set( float(value) )
         # fix vexp
         mej, ek = tailmodels.t0_to_Mej_Ek( self.var_t0.get(), self.var_vexp.get() )
         self.var_mej.set( mej )
         self.var_ek.set( ek )
         self.var_taum.set( arnettmodels.Mej_Ek_to_taum(self.var_mej.get(), self.var_ek.get()) )
         self.show_mbolmodel()
      l = Tk.Label(menu,text='Tail t0 (day) =',fg="black")      
      l.grid(row=6,column=0)
      self.var_t0 = Tk.DoubleVar(value=tailmodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()))
      spin = Tk.Spinbox(menu, textvariable=self.var_t0, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_t0, orient='horizontal',
                       length=200, command=update_t0)
      spin['to'] = constants.t0_bounds[1]
      spin['from'] = constants.t0_bounds[0]
      spin['increment'] = 0.1
      slide['to'] = constants.t0_bounds[1]
      slide['from'] = constants.t0_bounds[0]
      slide['digits'] = 4
      slide['resolution'] = 0.1
      spin.grid(row=6, column=1, sticky='news')
      slide.grid(row=6, column=2, sticky='news')
      
      # vexp
      def update_vexp(value):
         self.var_vexp.set( float(value) )
         # fix taum
         mej, ek = arnettmodels.taum_to_Mej_Ek( self.var_taum.get(), self.var_vexp.get() )
         self.var_mej.set( mej )
         self.var_ek.set( ek )
         self.var_t0.set( tailmodels.Mej_Ek_to_t0(self.var_mej.get(), self.var_ek.get()) )
         self.show_mbolmodel()
      l = Tk.Label(menu,text='V phot (10**3 km/s) =',fg="black")      
      l.grid(row=7,column=0)
      self.var_vexp = Tk.DoubleVar(value=arnettmodels.Mej_Ek_to_vej(self.var_mej.get(), self.var_ek.get()))
      spin = Tk.Spinbox(menu, textvariable=self.var_vexp, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_vexp, orient='horizontal',
                       length=200, command=update_vexp)
      spin['to'] = constants.vm_bounds[1]
      spin['from'] = constants.vm_bounds[0]
      spin['increment'] = 0.2
      slide['to'] = constants.vm_bounds[1]
      slide['from'] = constants.vm_bounds[0]
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
      self.var_me = Tk.DoubleVar(value=constants.me_p0)
      spin = Tk.Spinbox(menu, textvariable=self.var_me, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_me, orient='horizontal',
                       length=200, command=update_me)
      spin['to'] = constants.me_bounds[1]
      spin['from'] = constants.me_bounds[0]
      spin['increment'] = 0.01
      slide['to'] = constants.me_bounds[1]
      slide['from'] = constants.me_bounds[0]
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
      self.var_re = Tk.DoubleVar(value=constants.re_p0)
      spin = Tk.Spinbox(menu, textvariable=self.var_re, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_re, orient='horizontal',
                       length=200, command=update_re)
      spin['to'] = constants.re_bounds[1]
      spin['from'] = constants.re_bounds[0]
      spin['increment'] = 1
      slide['to'] = constants.re_bounds[1]
      slide['from'] = constants.re_bounds[0]
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
      self.var_ee = Tk.DoubleVar(value=np.log10(constants.ee_p0))
      spin = Tk.Spinbox(menu, textvariable=self.var_ee, wrap=True,
                        width=10, command=self.show_mbolmodel)
      slide = Tk.Scale(menu, variable=self.var_ee, orient='horizontal',
                       length=200, command=update_ee)      
      spin['to'] = np.log10(constants.ee_bounds[1])
      spin['from'] = np.log10(constants.ee_bounds[0])
      spin['increment'] = .1
      slide['to'] = np.log10(constants.ee_bounds[1])
      slide['from'] = np.log10(constants.ee_bounds[0])
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
         snerun.snobject_pars['multiband_main_xrangep'][0],
         snerun.snobject_pars['multiband_main_xrangep'][1],
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
                      flux/fm*snerun.snobject_pars['flux_scale'], color='r', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )         
      #if showarnett+showtail+showsbo>0:
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
      times_tail = np.arange(ts, max(constants.ts_bounds[1], 120), 1)
      self.Message('Mni=%.2f Mej=%.2f Ek=%.2f\n\nvexp=%.2f taum=%.1f t0=%.1f\n\ntexp=%.2f ts=%.2f\n\nMe=%.2f Re=%.1f Ee=%.3f'% (mni,mej,ek,vexp,taum,t0,texp,ts,me,re,ee))
      
      # plot
      if not 'wref' in self.__dict__: self.wref = []
      else:
         for w in self.wref:
            try: self.ax.lines.remove(w())
            except: pass
      if showarnett == 1:
         Larnett = arnettmodels.Arnett_fit_taum(times_arnett, mni, taum)         
         self.ax.plot(times_arnett+texp, Larnett, color='r', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )
      if showtail == 1:
         Ltail = tailmodels.tail_fit_t0(times_tail, mni, t0)
         self.ax.plot(times_tail+texp, Ltail, color='orange', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )
      if showsbo == 1:
         Lsbo = sbomodels.shock_fit(times_sbo, me, re, ee)
         self.ax.plot(times_sbo+texp, Lsbo, color='blue', ls='--')
         self.wref.append( ref(self.ax.lines[-1]) )
      if showarnett+showtail+showsbo>0:
         self.canvas.draw_idle()
   
   def Message(self,message,disable=True,restart=True):
      
      self.textbox.config(state='normal',foreground='black')
      if restart: self.textbox.delete(1.0,Tk.END)
      self.textbox.insert(Tk.END,message)
      self.msgwindow.config(command=self.textbox.yview)
      if disable: self.textbox.config(state='disabled')
      
   def quit(self):
      if messagebox.askokcancel("Quit", "Do you want to quit?"):
         self.parent.quit()
         self.parent.destroy()

class cMenu(Tk.Frame):
   
   def __init__(self, parent):
      
      Tk.Frame.__init__(self, parent)
      self.parent = parent

      frame = Tk.Frame(self.parent)   
      frame.grid(row=0,column=0,sticky=Tk.NW)

      self.file_path = Tk.StringVar()
      browsebutton = Tk.Button(frame, text="configure file", command=self.browse)
      browsebutton.grid(row=0, column=0,sticky=Tk.W) 
      
      quitBut = Tk.Button(frame,text="quit",fg="red", command=self.quit,)
      quitBut.grid(row=0,column=1,sticky=Tk.W)
      
   def browse(self):
      filetypes = [
         ('SDAPY configure file (*.default)', '*.txt'),         
      ]      
      filename = filedialog.askopenfilename(initialdir=os.getcwd(),
                                 title='file name', filetypes=filetypes)
      self.file_path.set(filename)
      
   def quit(self):
      if messagebox.askokcancel("Quit", "Do you want to quit?"): sys.exit()
      
################################   call to main #########################
def main():   
   root = Tk.Tk()
   screen_width = root.winfo_screenwidth()
   screen_height = root.winfo_screenheight()
   root.geometry("%dx%d+10+10"%(screen_width/1.5,screen_height/1.2))
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
