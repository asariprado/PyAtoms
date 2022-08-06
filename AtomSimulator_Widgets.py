 # -*- coding: utf-8 -*-
"""
ATOM SIMULATOR WIDGETS
Created on Mon Nov 15 14:45:06 2021
@author: Asari
"""

import sys
import os
import webbrowser

import numpy as np


from numpy import fft as npf
import math as m
from numpy import cos as cos
from numpy import sin as sin
from numpy import tan as tan
from numpy import arcsin as arcsin
from numpy import arccos as arccos
from numpy import arctan as arctan

from numpy import sqrt as sqrt
from numpy import log as log
from numpy import pi as pi



from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.image as mplimg

# PyQT imports for creating widgets/etc
from PyQt5.QtWidgets import *#QApplication, QCheckBox, QFileDialog, QMessageBox, QWidget, QLabel, QButtonGroup, QPushButton, QSpinBox, QMenu, QComboBox, QMainWindow, QHBoxLayout, QVBoxLayout, QSlider, QGroupBox, QGridLayout, QRadioButton, QDialog, QLineEdit, QInputDialog, QToolTip
from PyQt5.QtCore import *#Qt
from PyQt5.QtGui import *#QPainter, QColor, QAction


# THESE TWO ARE FOR EMBEDDING MATPLOTLIB PLOTS INTO PYQT5 GUIs
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# FigureCanvasQTAgg: It will provide the canvas for the figure
# NavigationToolbar2QT: It will provide the tool bar for the matplotlib figs (zooming in, panning, saving image, etc)
# https://www.geeksforgeeks.org/how-to-embed-matplotlib-graph-in-pyqt5/


from hexatoms import hexatoms
from squareatoms import squareatoms
from moirelattice import moirelattice


class SimulatorWidget(QWidget):

	def __init__(self, parent, x, y, w, h):
		super().__init__(parent)
		
		# Define all default values 
		self.moireBtn = "Single"
		self.pix = 256
		self.L = 7
	
		self.theta_im = 0
		self.theta_tw = 0
		self.theta_tw2 = 0 

		self.a = 0.3
		self.b = 0.3
		self.c = 0.3

		self.lattice1 = "Hexagonal"
		self.lattice2 = "Hexagonal"
		self.lattice3 = "Hexagonal"

		self.e11 = 0
		self.e12 = 0
		self.e22 = 0
		self.d11 = 0
		self.d12 = 0
		self.d22 = 0 
		self.f11 = 0
		self.f12 = 0
		self.f22 = 0 

		# self.honeycomb = 0
		# self.honeycomb2 = 0
		# self.honeycomb3 = 0
		self.origin1 = "Hollow"
		self.origin2 = "Hollow"
		self.origin3 = "Hollow"

		self.alpha1 = 1
		self.beta1 = 0
		self.alpha2 = 1
		self.beta2 = 0
		self.alpha3 = 1
		self.beta3 = 0

		self.filter_bool = False
		self.sigma = 0
		self.sigma_real = self.sigma*self.L/(self.pix-1)  *np.sqrt(2*np.log(2))

		self.saveFileName = ''

		self.colormap_RS = 'magma'
		self.colormap_FFT = 'Blues'
		self.colormapList = ['viridis', 'viridis_r', 'plasma', 'plasma_r', 'inferno', 'inferno_r', 'magma', 'magma_r','cividis', 'cividis_r',
							'Greys', 'Greys_r', 'Purples',  'Purples_r', 'Blues', 'Blues_r',  'Greens', 'Greens_r', 'Oranges', 'Oranges_r','Reds','Reds_r',
				            'YlOrBr','YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'OrRd', 'OrRd_r', 'PuRd', 'PuRd_r', 'RdPu', 'RdPu_r', 'BuPu', 'BuPu_r',
				            'GnBu', 'GnBu_r', 'PuBu', 'PuBu_r',  'YlGnBu', 'YlGnBu_r', 'PuBuGn', 'PuBuGn_r', 'BuGn', 'BuGn_r', 'YlGn', 'YlGn_r',
				            'binary','binary_r', 'gist_yarg','gist_yarg_r', 'gist_gray','gist_gray_r', 'gray','gray_r', 'bone','bone_r', 'pink','pink_r',
				            'spring', 'spring_r', 'summer', 'summer_r', 'autumn', 'autumn_r', 'winter', 'winter_r', 'cool', 'cool_r', 'Wistia', 'Wistia_r',
				            'hot', 'hot_r', 'afmhot', 'afmhot_r', 'gist_heat', 'gist_heat_r', 'copper', 'copper_r',
				            'PiYG',  'PiYG_r', 'PRGn','PRGn_r', 'BrBG', 'BrBG_r', 'PuOr', 'PuOr_r', 'RdGy','RdGy_r', 'RdBu','RdBu_r',
         				   	'RdYlBu', 'RdYlBu_r', 'RdYlGn','RdYlGn_r', 'Spectral', 'Spectral_r', 'coolwarm', 'coolwarm_r', 'bwr', 'bwr_r', 
         				   	'seismic','seismic_r',  'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r',  'hsv', 'hsv_r',
         				   	'Pastel1','Pastel1_r', 'Pastel2',  'Pastel2_r', 'Paired', 'Paired_r', 'Accent', 'Accent_r', 'Dark2', 'Dark2_r', 'Set1','Set1_r', 'Set2', 'Set2_r', 'Set3','Set3_r',
				            'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c','tab20c_r','flag', 'flag_r', 
				            'prism','prism_r', 'ocean', 'ocean_r', 'gist_earth', 'gist_earth_r', 'terrain', 'terrain_r', 'gist_stern', 'gist_stern_r',
				            'gnuplot', 'gnuplot_r', 'gnuplot2','gnuplot2_r', 'CMRmap', 'CMRmap_r', 'cubehelix', 'cubehelix_r', 'brg', 'brg_r',
				            'gist_rainbow', 'gist_rainbow_r', 'rainbow','rainbow_r', 'jet', 'jet_r', 'nipy_spectral', 'nipy_spectral_r','gist_ncar', 'gist_ncar_r']

		self.figure =plt.figure(figsize=(10,10))
		
		self.Z, self.fftZ = hexatoms(self.pix, self.L, self.a, self.theta_im, self.e11, self.e12, self.e22, self.alpha1, self.beta1, self.origin1)
		
		

		# STM topography time estimator
		self.vt = 20 # velocity of tip scanner, for calculating how long the image should take when doing actual stm measurements
		self.time_tot = (2 * self.pix * self.L) / self.vt # Total time in seconds
		self.hrs = int(self.time_tot // 3600)
		self.remain = self.time_tot - self.hrs*3600
		self.mins = int(self.remain // 60)
		self.sec = int(self.time_tot - (3600*self.hrs) - (60*self.mins))

		# dI/dV Map time estimator
		self.tps = 1 # time per spectra..
		self.time_map_tot = (self.tps*self.pix*self.pix) + self.time_tot #+ ((2* self.L) / (self.pix*self.vt))# Total time in seconds
		self.days_map = int(self.time_map_tot // 86400)
		self.remain_map = self.time_map_tot - (self.days_map*86400)
		self.hrs_map = int(self.remain_map // 3600)
		self.mins_map = int(self.time_map_tot - (self.days_map*86400) - (3600*self.hrs_map)) // 60
		

		# self.initMoireBtn()
		# self.initImageParameters()

		# self.initCalcWidget()
		# self.initMapCalcWidget()

		# self.initLattice1Parameters()
		# self.initLattice2Parameters()
		# self.initLattice3Parameters()

		# self.initColormapDropdown()
		# self.initMatplotlibFig()
	
		# self.plotAtoms() # Uncommenting this makes an error because of the self.figure.clear() part...

	def initImageParameters(self):
		groupBox = QGroupBox("Image parameters")
		vlayout = QVBoxLayout(self)

		# Pixels
		self.pix_label = QLabel("Pixels:")
		self.pix_label.setToolTip("Create a square (pix x pix) image")
		self.pix_input = QLineEdit(self)
		self.pix_input.returnPressed.connect(self.updatePix) # Connect this intput dialog whenever the enter/return button is pressed
		## !!!! use # .editingFinished.connect instead of .returnPressed.connect because it detects when enter/tab is pressed or if you click on a different widget in the gui. then it'll update  
		## From https://doc.qt.io/qt-6/qlineedit.html#editingFinished
		# self.pix_input.editingFinished.connect(self.updatePix) # nvm it doesnt work correctly
		self.pix_input.setPlaceholderText(str(self.pix))
		self.pix_input.setFixedWidth(100)
		self.pix_input.setToolTip("Create a square (pix x pix) image")

		self.pix_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update!
		self.pix_btn.clicked.connect(self.updatePix)
		self.pix_btn.setAutoDefault(False)

		# Define label to display the value of the slider next to the textbox
		self.pix_SliderLabel = QLabel(" pix", self)
		# self.pix_SliderLabel.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.pix_SliderLabel.setMinimumWidth(30)

		hbox = QHBoxLayout(self)
		hbox.addWidget(self.pix_label)
		hbox.addWidget(self.pix_input)
		hbox.addWidget(self.pix_SliderLabel)
		hbox.addWidget(self.pix_btn)

		vlayout.addLayout(hbox)
		vlayout.setSpacing(3) # To reduce the spacing (vertical?) between the widgets. saves so much space


		# Length of image
		self.length_Label = QLabel("L: ")
		self.length_Label.setToolTip("Length of image in nanometers") 
		self.L_input = QLineEdit(self)
		self.L_input.returnPressed.connect(self.updateL) # Connect this intput dialog whenever the enter/return button is pressed
		# self.L_input.returnPressed.connect(self.update_vt) # Connect to update_vt function to calculate the new estimate time of STM image bc it depends on the length of image
		self.L_input.setPlaceholderText(str(self.L))
		self.L_input.setFixedWidth(100)
		self.L_input.setToolTip("Length of image in nanometers") 

		# Literally just creates the little label next to the text box that says the units
		self.L_Label = QLabel(" nm", self) # Display the corrected value. only up to 2 decimal pts
		# self.L_Label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.L_Label.setMinimumWidth(30)

		self.L_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.L_btn.clicked.connect(self.updateL)
		self.L_btn.setAutoDefault(False)

		hbox2 = QHBoxLayout(self)
		hbox2.addWidget(self.length_Label)
		hbox2.addWidget(self.L_input)
		hbox2.addWidget(self.L_Label)
		hbox2.addWidget(self.L_btn)

		vlayout.addLayout(hbox2)



		# Image twist angle theta
		self.theta_label = QLabel("Theta:")
		self.theta_label.setToolTip("Rotates the entire image by the input offset angle") 
		self.theta_im_input = QLineEdit(self)
		self.theta_im_input.returnPressed.connect(self.updateTheta) # Connect this intput dialog whenever the enter/return button is pressed
		self.theta_im_input.setPlaceholderText(str(self.theta_im))
		self.theta_im_input.setFixedWidth(100)
		self.theta_im_input.setToolTip("Rotates the entire image by the input offset angle")

		# Literally just creates the little label next to the text box that says the units
		self.thetaLabel = QLabel('deg', self)
		# self.thetaLabel.setAlignment(Qt.AlignCenter | Qt.AlignVCenter) # Idk what this does
		self.thetaLabel.setMinimumWidth(30)

		self.theta_btn = QPushButton("Rotate", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.theta_btn.clicked.connect(self.updateTheta)
		self.theta_btn.setAutoDefault(False)

		# Place the text box and display value text next to each other horizontally
		hbox3 = QHBoxLayout()
		hbox3.addWidget(self.theta_label)
		hbox3.addWidget(self.theta_im_input)
		# hbox3.addSpacing(10)
		hbox3.addWidget(self.thetaLabel)
		hbox3.addWidget(self.theta_btn)

		vlayout.addLayout(hbox3) # Add the widgets and value text to the groupbox

		groupBox.setLayout(vlayout)

		return groupBox

	def updatePix(self):

		# if self.pix_input.hasFocus(): # to only try updating if you are clicked on/typing in the pix_input box. bc otherwise, if you put an invalid input, error messages keep popping up even after you click away from it.. its annoying.. jk it doesnt work :( onnly if you press enter it works, not if you do tab / click away)
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.pix = eval(self.pix_input.text()) #float(self.a_input.text())
			
			if self.pix >= 8192:
				raise ValueError
			else: 
				self.pix_SliderLabel.setText('pix') 
				self.pix_input.setPlaceholderText(str(self.pix))
				self.update_calc()
				self.update_map_calc()
				self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.pix_error = QMessageBox()
			self.pix_error.setWindowTitle("Error")
			# if self.pix >= 8192:
			# self.pix_error.setText('Number of pixels exceeds 8192. Please enter a smaller number.')
			# else:
			self.pix_error.setText("Type in a valid number or numerical expression. Make sure the number of pixels does not exceed 8192")
			
			self.pix_error.setInformativeText("Your input for pix is: " +  str(self.pix_input.text()))
			self.pix_error.setIcon(QMessageBox.Warning)
			self.pix_error.setStandardButtons(QMessageBox.Retry)
			x = self.pix_error.exec()

	def updateL(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.L = eval(self.L_input.text()) #float(self.a_input.text())
			self.L_input.setPlaceholderText(str(self.L))
			self.L_Label.setText('nm') 
			self.update_calc()
			self.update_map_calc()
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.L_error = QMessageBox()
			self.L_error.setWindowTitle("Error")
			self.L_error.setText("Type in a number or numerical expression.")
			self.L_error.setInformativeText("Your input for L is: " +  str(self.L_input.text()))
			self.L_error.setIcon(QMessageBox.Warning)
			self.L_error.setStandardButtons(QMessageBox.Retry)
			x = self.L_error.exec()

	def updateTheta(self):
		try: 
			self.theta_im = eval(self.theta_im_input.text())
			self.thetaLabel.adjustSize()
			self.theta_im_input.setPlaceholderText(str(self.theta_im))
			self.plotAtoms() 
		except:
			self.theta_im_error = QMessageBox()
			self.theta_im_error.setWindowTitle("Error")
			self.theta_im_error.setText("Type in a number or numerical expression")
			self.theta_im_error.setInformativeText("Your input for theta_im is: " +  str(self.theta_im_input.text()))
			self.theta_im_error.setIcon(QMessageBox.Warning)
			self.theta_im_error.setStandardButtons(QMessageBox.Retry)
			x = self.theta_im_error.exec()

	def initMoireBtn(self):
		groupBox = QGroupBox("Moire lattice?")
		groupBox.setToolTip("Choose whether to create a moire pattern with two sublattices, or plot only a single lattice")

		hbox = QHBoxLayout()

		self.noMoire = QRadioButton("Single")
		self.noMoire.setToolTip("Create a single lattice using Lattice 1 parameters")
		self.yesMoire = QRadioButton("Bilayer")
		self.yesMoire.setToolTip("Creates a moire pattern with two overlayed lattices")
		self.trilayer = QRadioButton("Trilayer")
		self.trilayer.setToolTip("Create a trilayer moire lattice")

		self.noMoire.setChecked(True) # Default value is no moire lattice. plot only 1 lattice
		

		## IMPORTANTTTT: have to connect all buttons to the SAME update__button function so theyll be like mutually exclusive. 
		## so if you click one, it means the others are set to false, etc. this is defined in the updating moire btn function
		self.yesMoire.toggled.connect(self.updateMoireBtn)
		self.noMoire.toggled.connect(self.updateMoireBtn)
		self.trilayer.toggled.connect(self.updateMoireBtn)



		# Create QButtonGroup to be mutually exclusive  buttons, from https://www.programcreek.com/python/example/108083/PyQt5.QtWidgets.QButtonGroup
		self.moire_btn_group = QButtonGroup(self)
		self.moire_btn_group.addButton(self.noMoire)
		self.moire_btn_group.addButton(self.yesMoire)
		self.moire_btn_group.addButton(self.trilayer)


		hbox.addWidget(self.noMoire)
		hbox.addWidget(self.yesMoire)
		hbox.addWidget(self.trilayer)

		groupBox.setLayout(hbox)
		groupBox.setContentsMargins(0,0,0,0) # Sets the left , top , right , and bottom margins to use around the layout.


		return groupBox

	def updateMoireBtn(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == "Bilayer":
				self.moireBtn = "Bilayer"
			elif radio_btn.text() == "Single":
				self.moireBtn = "Single"
			elif radio_btn.text() == "Trilayer":
				self.moireBtn = "Trilayer"

			# Update the plot every time user changes the button choice
		self.plotAtoms() 

	def initCalcWidget(self):
		groupBox = QGroupBox("SPM Time Estimator")
		groupBox.setToolTip("Calculate estimated time to take STM topography")
		vlayout = QVBoxLayout(self)

		self.calc_topo_btn = QCheckBox("Calculate     ")
		self.calc_topo_btn.setChecked(False)
		self.calc_topo_btn.stateChanged.connect(self.update_calc)
		self.calc_topo_btn.setToolTip("Check to calculate")
		# self.calc_btn_grp = QButtonGroup()
		# self.calc_btn_grp.addButton(self.calc_topo_btn, 1)

		hbox = QHBoxLayout(self)
		hbox.addWidget(self.calc_topo_btn)

		vlayout = QVBoxLayout(self)
		# vlayout.addLayout(hbox)

		# Create vt input text box
		self.vt_input = QLineEdit(self)
		self.vt_input.returnPressed.connect(self.update_vt) # Connect this intput dialog whenever the enter/return button is pressed
		# self.vt_input.returnPressed.connect(self.update_tps)
		self.vt_input.setPlaceholderText(str(self.vt))
		self.vt_input.setFixedWidth(60)

		self.vt_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update!
		self.vt_btn.clicked.connect(self.update_vt)
		self.vt_btn.setAutoDefault(False)

		# Define label to display the value of the slider next to the textbox
		self.vt_label = QLabel("Tip speed: ", self)
		self.vt_label.setToolTip("Velocity of tip")
		# self.vt_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.vt_label.setMinimumWidth(30)
		self.vt_nm_label = QLabel(" nm/s", self)
		# self.vt_nm_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.vt_nm_label.setMinimumWidth(30)


		self.topo_eta_label = QLabel("Est. time: " + str(self.hrs) + 'h ' + str(self.mins) + 'min ' + str(self.sec) + 's')

		hbox.addWidget(self.topo_eta_label)


		hbox2 = QHBoxLayout(self)
		hbox2.addWidget(self.vt_label)
		hbox2.addWidget(self.vt_input)
		hbox2.addWidget(self.vt_nm_label)
		hbox2.addWidget(self.vt_btn)

		vlayout.addLayout(hbox)
		vlayout.addLayout(hbox2)

		vlayout.setSpacing(1)

		groupBox.setLayout(vlayout)

		return groupBox

	def update_calc(self):
		if (self.calc_topo_btn).isChecked():
			# Calculate the estimated time it will take to take an STM topograph image
			# with the given pix, L, velocity of tip scanner vt
			self.time_tot = (2 * self.pix * self.L) / self.vt # Total time in seconds

			self.hrs = int(self.time_tot // 3600)
			self.remain = self.time_tot - self.hrs*3600
			self.mins = int(self.remain // 60)
			self.sec = int(self.time_tot - (3600*self.hrs) - (60*self.mins))

			self.topo_eta_label.setText("Est. time: " + str(self.hrs) + 'h ' + str(self.mins) + 'min ' + str(self.sec) + 's')
		else:
			pass
		
	def update_vt(self):

		# if self.vt_input.hasFocus(): # commenting out bc it doesnt work like i wanted it to :(

		if self.calc_topo_btn.isChecked():
			try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
					# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
				self.vt = eval(self.vt_input.text()) 
				self.vt_input.setPlaceholderText(str(self.vt))

				self.update_calc()
				self.update_map_calc()

			except:
				# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
				# https://www.w3schools.com/python/python_try_except.asp
				# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
				self.vt_error = QMessageBox()
				self.vt_error.setWindowTitle("Error")
				self.vt_error.setText("Type in a number or numerical expression.")
				self.vt_error.setInformativeText("Your input for vt is: " +  str(self.vt_input.text()))
				self.vt_error.setIcon(QMessageBox.Warning)
				self.vt_error.setStandardButtons(QMessageBox.Retry)
				x = self.vt_error.exec()
		else:
			pass

	def initMapCalcWidget(self):
		groupBox = QGroupBox("Spectroscopy Map Time Estimator")
		groupBox.setToolTip("Calculate estimated time to take dI/dV map")
		vlayout = QVBoxLayout(self)

		self.calc_map_btn = QCheckBox("Calculate     ")
		self.calc_map_btn.setChecked(False)
		self.calc_map_btn.stateChanged.connect(self.update_map_calc)
		self.calc_map_btn.setToolTip("Check to calculate")
		# self.calc_btn_grp = QButtonGroup()
		# self.calc_btn_grp.addButton(self.calc_topo_btn, 1)

		hbox = QHBoxLayout(self)
		hbox.addWidget(self.calc_map_btn)

		vlayout = QVBoxLayout(self)
		# vlayout.addLayout(hbox)

		# Create tps (time per spectra) input text box
		self.tps_input = QLineEdit(self)
		self.tps_input.returnPressed.connect(self.update_tps) # Connect this intput dialog whenever the enter/return button is pressed
		self.tps_input.setPlaceholderText(str(self.tps))
		self.tps_input.setFixedWidth(40)

		self.tps_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update!
		self.tps_btn.clicked.connect(self.update_tps)
		self.tps_btn.setAutoDefault(False)

		# Define label to display the value of the slider next to the textbox
		self.tps_label = QLabel("Time per spectra: ", self)
		self.tps_label.setToolTip("Time per spectra")
		# self.tps_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.tps_label.setMinimumWidth(30)
		self.tps_s_label = QLabel(" s", self)
		# self.tps_s_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.tps_s_label.setMinimumWidth(30)

		self.map_eta_label = QLabel("Est. time: " + str(self.days_map) + 'd' + str(self.hrs_map) + 'h' + str(self.mins_map) + 'm')

		# hbox3 = QHBoxLayout(self)
		hbox.addWidget(self.map_eta_label)


		hbox2 = QHBoxLayout(self)
		hbox2.addWidget(self.tps_label)
		hbox2.addWidget(self.tps_input)
		hbox2.addWidget(self.tps_s_label)
		hbox2.addWidget(self.tps_btn)

		vlayout.addLayout(hbox)
		vlayout.addLayout(hbox2)

		vlayout.setSpacing(3)

		groupBox.setLayout(vlayout)

		return groupBox

	def update_map_calc(self):
		if (self.calc_map_btn).isChecked():
			# Calculate the estimated time it will take to take a dI/dV map 
			self.time_map_tot = (self.tps*self.pix*self.pix) + self.time_tot # This is the total time in seconds

			# Convert time in seconds to days, hours, mins:
			self.days_map = int(self.time_map_tot // (60*60*24)) # Calculate how many full days in the total time: mod by 60*60*24 bc 60 seconds in a minute, 60 mins in an hour, 24 hours in a day
			self.remain_map = self.time_map_tot - (self.days_map*86400) # Calculate time remaining after subtracting the hours we already accounted for
			self.hrs_map = int(self.remain_map // (60*60)) # Calculate the hours by modding by 60*60, bc 60 secs in a min, 60 min in an hour
			self.mins_map = int(self.time_map_tot - (self.days_map*60*60*24) - (60*60*self.hrs_map)) // 60 # Calculate the number of minutes after subtracting for the days/hours we already accounted for. mod by 60  bc 60 secs in a minute,
			## These calculations are correct, verified I got the right seconds --> days/hrs/mins from this website converter https://www.satsig.net/training/seconds-days-hours-minutes-calculator.htm 
			
			# self.hrs_map = int(self.time_map_tot // (60*60)) # converting to hours . 3600 sec in an hour
			# self.remain_map = self.time_map_tot - self.hrs_map*3600
			# self.mins_map = int(self.remain_map // 60)
			# self.sec_map = int(self.time_map_tot - (3600*self.hrs_map) - (60*self.mins_map))
			# self.topo_calc_params.setText("pix: " + str(self.pix) + '\nL: ' + str(self.L) +'\nvt: ' + str(self.vt))
			# self.map_eta_label.setText("Est. time: " + str(self.hrs_map) + 'h ' + str(self.mins_map) + 'min ' + str(self.sec_map) + 's')

			self.map_eta_label.setText("Est. time: " + str(self.days_map) + 'd' + str(self.hrs_map) + 'h' + str(self.mins_map) + 'm')

		else:
			pass
		
	def update_tps(self):

		# if self.vt_input.hasFocus(): # commenting out bc it doesnt work like i wanted it to :(

		if self.calc_map_btn.isChecked():
			try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
					# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
				self.tps = eval(self.tps_input.text()) #float(self.a_input.text())
				self.tps_input.setPlaceholderText(str(self.tps))

				self.update_map_calc()
				# total time estimated in hrs.mins.sec format
				# print(str(hrs) + 'h ' + str(mins) + 'min ' + str(sec) + 's')

			except:
				# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
				# https://www.w3schools.com/python/python_try_except.asp
				# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
				self.tps_error = QMessageBox()
				self.tps_error.setWindowTitle("Error")
				self.tps_error.setText("Type in a number or numerical expression.")
				self.tps_error.setInformativeText("Your input for time per spectra is: " +  str(self.tps_input.text()))
				self.tps_error.setIcon(QMessageBox.Warning)
				self.tps_error.setStandardButtons(QMessageBox.Retry)
				x = self.tps_error.exec()
		else:
			pass

	def initLattice1Parameters(self):

		lat1groupBox = QGroupBox("Lattice 1")


		# Add tabs to save space. 1 tab for lattice parameters (symmetry, periodicity, honeycomb/dots), the other for strain
		# https://pythonspot.com/pyqt5-tabs/
		self.lat1tabs = QTabWidget(self)

		self.tab1a = QWidget(self)
		self.tab2a = QWidget(self)
		self.tab3a = QWidget(self)
		# self.tab4a = QWidget(self)
		self.lat1tabs.addTab(self.tab1a, "Parameters")
		self.lat1tabs.addTab(self.tab3a, "Sublattices")
		# self.lat1tabs.addTab(self.tab4a, "Origin")
		self.lat1tabs.addTab(self.tab2a, "Strain")

		self.tab3a.setToolTip("Only works for hexagonal lattices")
		



		# # # # # # # # # # # # # # # # # # # # 
		### Pick lattice 1 symmetry ###
		self.symm_label = QLabel("Symmetry:", self)
		self.symm_label.setToolTip("Choose symmetry of first lattice")
		self.hex1btn = QRadioButton("Hexagonal")
		self.sq1btn = QRadioButton("Square")
		self.hex1btn.setChecked(True)
		self.hex1btn.setToolTip("Create hexagonal lattice")
		self.sq1btn.setToolTip("Create square lattice")

		# Connect btn to update functions when clicked https://www.tutorialspoint.com/pyqt/pyqt_qpushbutton_widget.htm
		self.hex1btn.toggled.connect(self.updateLattice1)
		self.sq1btn.toggled.connect(self.updateLattice1)

		# Create QButtonGroup to be mutually exclusive w the honeycomb buttons, from https://www.programcreek.com/python/example/108083/PyQt5.QtWidgets.QButtonGroup
		self.symm_btn_group = QButtonGroup(self)
		self.symm_btn_group.addButton(self.hex1btn)
		self.symm_btn_group.addButton(self.sq1btn)

		h1box = QHBoxLayout(self)
		h1box.addWidget(self.hex1btn)
		h1box.addWidget(self.sq1btn)

		# Add symmetry button widgets
		self.tab1a.layout = QVBoxLayout(self)
		self.tab1a.layout.addWidget(self.symm_label)
		self.tab1a.layout.addLayout(h1box)


		# # # # # # # # # # # # # # # # # # # # 
		### Pick lattice 1 periodicity a ###
		self.a_param_label = QLabel("Lattice constant (nm)", self)
		self.a_param_label.setToolTip("Periodicity of 1st lattice (spacing between atoms) in nm")
		self.a_label = QLabel("a:", self)
		self.a_label.setToolTip("Periodicity of 1st lattice (spacing between atoms) in nm")
		self.a_input = QLineEdit(self)
		self.a_input.returnPressed.connect(self.update_a) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.a_input.setPlaceholderText(str(self.a))
		self.a_input.setToolTip("Periodicity of 1st lattice (spacing between atoms) in nm")
		self.a_input.setMinimumWidth(85)

		# Define label to display the value of the slider next to the slider
		self.a_nm_Label = QLabel(" nm", self) # Display the corrected value. only up to 2 decimal pts
		# self.a_nm_Label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.a_nm_Label.setMinimumWidth(30)

		self.a_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.a_btn.clicked.connect(self.update_a)
		self.a_btn.setAutoDefault(False)


		h2box = QHBoxLayout(self)
		h2box.addWidget(self.a_label)
		h2box.addWidget(self.a_input)
		h2box.addWidget(self.a_nm_Label)
		h2box.addWidget(self.a_btn)


		# Add lattice constant a widgets to the tab 1
		self.tab1a.layout.addWidget(self.a_param_label)
		self.tab1a.layout.addLayout(h2box)


		



		# # # # # # # # # # # # # # # # # # # # 
		### Pick strain on lattice 1 ###
		self.strain1_label = QLabel("Strain (in percent, i.e. 3%)")

		# e11 input
		self.e11_input = QLineEdit(self)
		self.e11_input.setPlaceholderText(str(self.e11))
		self.e11_input.setFixedWidth(60)
		self.e11_label = QLabel("e11", self)
		self.e11_SliderLabel = QLabel(str("%"), self) 
		self.e11_input.returnPressed.connect(self.update_e11)
		self.e11_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.e11_btn.clicked.connect(self.update_e11)
		self.e11_btn.setAutoDefault(False)
		self.e11_label.setToolTip("Strain tensor elements")
		self.e11_input.setToolTip("Strain tensor elements")

		# e12 input
		self.e12_input = QLineEdit(self)
		self.e12_input.setPlaceholderText(str(self.e12))
		self.e12_input.setFixedWidth(60)
		self.e12_label = QLabel("e12", self)
		self.e12_SliderLabel = QLabel(str("%"), self) 
		self.e12_input.returnPressed.connect(self.update_e12)
		self.e12_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.e12_btn.clicked.connect(self.update_e12)
		self.e12_btn.setAutoDefault(False)
		self.e12_label.setToolTip("Strain tensor elements")
		self.e12_input.setToolTip("Strain tensor elements")

		# e22 input
		self.e22_input = QLineEdit(self)
		self.e22_input.setPlaceholderText(str(self.e22))
		self.e22_input.setFixedWidth(60)
		self.e22_label = QLabel("e22", self)
		self.e22_SliderLabel = QLabel(str("%"), self) 
		self.e22_input.returnPressed.connect(self.update_e22)
		self.e22_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.e22_btn.clicked.connect(self.update_e22)
		self.e22_btn.setAutoDefault(False)
		self.e22_label.setToolTip("Strain tensor elements")
		self.e22_input.setToolTip("Strain tensor elements")


		h3box = QHBoxLayout(self)
		h3box.addWidget(self.e11_label)
		h3box.addWidget(self.e11_input)
		h3box.addWidget(self.e11_SliderLabel)
		h3box.addWidget(self.e11_btn)


		h4box = QHBoxLayout(self)
		h4box.addWidget(self.e12_label)
		h4box.addWidget(self.e12_input)
		h4box.addWidget(self.e12_SliderLabel)
		h4box.addWidget(self.e12_btn)
		

		h5box = QHBoxLayout(self)
		h5box.addWidget(self.e22_label)
		h5box.addWidget(self.e22_input)
		h5box.addWidget(self.e22_SliderLabel)
		h5box.addWidget(self.e22_btn)


		# Put all the strain widgets into tab2 
		self.tab2a.layout = QVBoxLayout(self)
		self.tab2a.layout.addWidget(self.strain1_label)
		self.tab2a.layout.addLayout(h3box)
		self.tab2a.layout.addLayout(h4box)
		self.tab2a.layout.addLayout(h5box)





		# # # # # # # # # # # # # # # # # # # # # 
		# ### Pick origin site for lattice 1 ###
		self.origin1label = QLabel("Lattice site at origin")
		
		# Create radiobuttons
		self.hollowsite1 = QRadioButton("Hollow") 
		self.Asite1 = QRadioButton("A-site")
		self.Bsite1 = QRadioButton("B-site")

		self.hollowsite1.setChecked(True)

		# Set hover tool tips
		self.origin1label.setToolTip("Only works for hexagonal lattices.")
		self.hollowsite1.setToolTip("No atom at the origin")
		self.Asite1.setToolTip("A-atom at the origin")
		self.Bsite1.setToolTip("B-atom at the origin")

		# # Create QButtonGroup to be mutually exclusive  buttons. This is all u have to do for them to work correctly! 
		self.origin1_group = QButtonGroup(self)
		self.origin1_group.addButton(self.hollowsite1)
		self.origin1_group.addButton(self.Asite1)
		self.origin1_group.addButton(self.Bsite1)

		# ### IMPORTANTTTT: have to connect all buttons to the same update__button function so theyll be like mutually exclusive. so if you click one, it means the others are set to false, etc. this is defined in the updating function
		self.hollowsite1.toggled.connect(self.updateOrigin1)
		self.Asite1.toggled.connect(self.updateOrigin1)
		self.Bsite1.toggled.connect(self.updateOrigin1)

		h6box = QHBoxLayout(self)
		h6box.addWidget(self.hollowsite1)
		h6box.addWidget(self.Asite1)
		h6box.addWidget(self.Bsite1)






		self.sublattices_label1 = QLabel("Pick strength of sublattices Z = A + B")
		self.sublattices_label2 = QLabel("Honeycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")


		self.sublattices_label = QLabel("\nWeight of sublattices")
		self.sublattices_label.setToolTip("Only works if hexagonal symmetry is selected.\nHoneycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")



		# # # # # # # # # # # # # # # # # # # # 
		### Pick alpha1 weight of sublattice a ###
		# self.alpha1_param_label = QLabel("Lattice constant (nm)", self)
		# self.alpha1_param_label.setToolTip("Weight of sublattice a")
		self.alpha1_label = QLabel("alpha1:", self)
		self.alpha1_label.setToolTip("Weight of sublattice A")
		self.alpha1_input = QLineEdit(self)
		self.alpha1_input.returnPressed.connect(self.update_alpha1) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.alpha1_input.setPlaceholderText(str(self.alpha1))
		self.alpha1_input.setToolTip("Weight of sublattice A")
		self.alpha1_input.setFixedWidth(60)

		
		self.alpha1_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.alpha1_btn.clicked.connect(self.update_alpha1)
		self.alpha1_btn.setAutoDefault(False)
		self.alpha1_btn.setToolTip("Weight of sublattice A")


		h7box = QHBoxLayout(self)
		h7box.addWidget(self.alpha1_label)
		h7box.addWidget(self.alpha1_input)
		h7box.addWidget(self.alpha1_btn)



		# # # # # # # # # # # # # # # # # # # # 
		### Pick beta1 weight of sublattice b ###
		# self.beta1_param_label = QLabel("Lattice constant (nm)", self)
		# self.beta1_param_label.setToolTip("Weight of sublattice B")
		self.beta1_label = QLabel("beta1:", self)
		self.beta1_label.setToolTip("Weight of sublattice B")
		self.beta1_input = QLineEdit(self)
		self.beta1_input.returnPressed.connect(self.update_beta1) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.beta1_input.setPlaceholderText(str(self.beta1))
		self.beta1_input.setToolTip("Weight of sublattice B")
		self.beta1_input.setFixedWidth(60)

		
		self.beta1_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.beta1_btn.clicked.connect(self.update_beta1)
		self.beta1_btn.setAutoDefault(False)
		self.beta1_btn.setToolTip("Weight of sublattice B")


		h8box = QHBoxLayout(self)
		h8box.addWidget(self.beta1_label)
		h8box.addWidget(self.beta1_input)
		h8box.addWidget(self.beta1_btn)






		

		# Add honeycomb/dots button widgets to the tab 3
		self.tab3a.layout = QVBoxLayout(self)
		self.tab3a.layout.addWidget(self.origin1label)
		self.tab3a.layout.addLayout(h6box)
		# self.tab3a.layout.addWidget(self.sublattices_label1)
		# self.tab3a.layout.addWidget(self.sublattices_label2)
		self.tab3a.layout.addWidget(self.sublattices_label)
		self.tab3a.layout.addLayout(h7box)
		self.tab3a.layout.addLayout(h8box)

		

		

		# Set the layouts for each tab...
		self.tab1a.setLayout(self.tab1a.layout)
		self.tab1a.layout.setSpacing(1)
		self.tab2a.setLayout(self.tab2a.layout)
		self.tab2a.layout.setSpacing(1)
		self.tab3a.setLayout(self.tab3a.layout)
		self.tab3a.layout.setSpacing(1)


		vlayout = QVBoxLayout(self)
		vlayout.addWidget(self.lat1tabs)
		lat1groupBox.setLayout(vlayout)
		vlayout.setSpacing(1)

		# lat1groupBox.setContentsMargins(0,15,0,0) # Sets the left , top , right , and bottom margins to use around the layout.


		return lat1groupBox

	def initLattice2Parameters(self):
		lat2groupBox = QGroupBox("Lattice 2")
		# lat2groupbox.setToolTip("Parameters for 2nd lattice when simulating a Moire lattice.\nChanging these values will not affect the image unless the 'Yes' button in Moire lattice is selected.")

		# Add tabs to save space. 1 tab for lattice parameters (symmetry, periodicity, honeycomb/dots), the other for strain
		# https://pythonspot.com/pyqt5-tabs/
		self.lat2tabs = QTabWidget(self)

		self.tab1b = QWidget(self)
		self.tab2b = QWidget(self)
		self.tab3b = QWidget(self)
		self.lat2tabs.addTab(self.tab1b, "Parameters")
		self.lat2tabs.addTab(self.tab3b, "Sublattices")
		self.lat2tabs.addTab(self.tab2b, "Strain")

			
		# # # # # # # # # # # # # # # # # # # # 
		### Pick lattice 2 symmetry ###
		self.symm_label = QLabel("Symmetry:", self)
		self.symm_label.setToolTip("Choose symmetry of second lattice")
		self.hex2btn = QRadioButton("Hexagonal")
		self.sq2btn = QRadioButton("Square")
		self.hex2btn.setChecked(True)
		self.hex2btn.setToolTip("Create hexagonal lattice")
		self.sq2btn.setToolTip("Create square lattice")

		# Connect btn to update functions when clicked https://www.tutorialspoint.com/pyqt/pyqt_qpushbutton_widget.htm
		self.hex2btn.toggled.connect(self.updateLattice2)
		self.sq2btn.toggled.connect(self.updateLattice2)

		# Create QButtonGroup to be mutually exclusive w the honeycomb buttons, from https://www.programcreek.com/python/example/108083/PyQt5.QtWidgets.QButtonGroup
		self.symm2_btn_group = QButtonGroup(self)
		self.symm2_btn_group.addButton(self.hex2btn)
		self.symm2_btn_group.addButton(self.sq2btn)

		h1box = QHBoxLayout(self)
		h1box.addWidget(self.hex2btn)
		h1box.addWidget(self.sq2btn)

		# Add symmetry button widgets to tab
		self.tab1b.layout = QVBoxLayout(self)
		self.tab1b.layout.addWidget(self.symm_label)
		self.tab1b.layout.addLayout(h1box)


		# # # # # # # # # # # # # # # # # # # # 
		### Pick lattice 2 periodicity b ###
		self.b_param_label = QLabel("\nLattice constant (nm)", self)
		self.b_label = QLabel("b:", self)
		# self.b_label.setFixedWidth(20)
		self.b_input = QLineEdit(self)
		self.b_input.returnPressed.connect(self.update_b) # Connect this intput dialog whenever the enter/return button is pressed
		self.b_input.setPlaceholderText(str(self.b))
		self.b_param_label.setToolTip("Periodicity of 2nd lattice (spacing between atoms) in nm")
		self.b_label.setToolTip("Periodicity of 2nd lattice (spacing between atoms) in nm")
		self.b_input.setToolTip("Periodicity of 2nd lattice (spacing between atoms) in nm")
		# self.b_label.setBuddy(self.b_input)
		self.b_input.setMinimumWidth(85)

		self.b_label.setMinimumHeight(15)
		self.b_input.setMinimumHeight(15)
		self.b_param_label.setMinimumHeight(15)

		# Define label to display the value of the slider next to the slider
		self.b_nm_Label = QLabel(" nm", self) # Display the corrected value. only up to 2 decimal pts
		# self.b_nm_Label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.b_nm_Label.setMinimumWidth(30)
		self.b_nm_Label.setMinimumHeight(15)

		self.b_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.b_btn.clicked.connect(self.update_b)
		self.b_btn.setAutoDefault(False)
		self.b_btn.setFixedWidth(60)
		self.b_btn.setMinimumHeight(15)

		h2box = QHBoxLayout(self)
		h2box.addWidget(self.b_label)
		h2box.addWidget(self.b_input)
		h2box.addWidget(self.b_nm_Label)
		h2box.addWidget(self.b_btn)

		# Add lattice constant a widgets to the tab 1
		self.tab1b.layout.addWidget(self.b_param_label)
		self.tab1b.layout.addLayout(h2box)



		# # # # # # # # # # # # # # # # # # # # 
		### Twist angle ###
		self.thetatw_title = QLabel('Twist angle', self)
		self.thetatw_input = QLineEdit(self)
		self.thetatw_input.returnPressed.connect(self.updateThetaTw) # Connect this intput dialog whenever the enter/return button is pressed
		self.thetatw_input.setPlaceholderText(str(self.theta_tw))
		self.thetatw_input.setToolTip("Relative twist angle between 1st & 2nd lattice")
		self.thetatw_input.setFixedWidth(60)
		self.thetatw_input.setMinimumHeight(15)

		# Define label to display the value of the slider next to the slider
		self.thetatwLabel = QLabel(' deg', self)
		self.thetatwLabel.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		# self.thetatwLabel.setMinimumWidth(30)
		self.thetatwLabel.setFixedWidth(30)
		self.thetatwLabel.setMinimumHeight(15)

		self.theta_tw_btn = QPushButton("Twist", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.theta_tw_btn.clicked.connect(self.updateThetaTw)
		self.theta_tw_btn.setAutoDefault(False)
		self.theta_tw_btn.setMinimumHeight(15)

		h7box = QHBoxLayout(self)
		h7box.addWidget(self.thetatw_input)
		h7box.addWidget(self.thetatwLabel)
		h7box.addWidget(self.theta_tw_btn)


		# Add twist angle widgets to the tab 1
		self.tab1b.layout.addWidget(self.thetatw_title)
		self.tab1b.layout.addLayout(h7box)



		# # # # # # # # # # # # # # # # # # # # 
		### Pick strain on lattice 2 ###
		self.strain2_label = QLabel("Strain (in percent, i.e. 3%)")

		# d11 input
		self.d11_input = QLineEdit(self)
		self.d11_input.setPlaceholderText(str(self.d11))
		self.d11_input.setFixedWidth(60) # Shorten the width of the text box
		self.d11_label = QLabel("d11", self)
		# self.d11_label.setFixedWidth(20)
		self.d11_SliderLabel = QLabel(str("%"), self) 
		# self.d11_SliderLabel.setFixedWidth(20)
		self.d11_input.returnPressed.connect(self.update_d11)
		self.d11_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.d11_btn.clicked.connect(self.update_d11)
		self.d11_btn.setAutoDefault(False)
		self.d11_label.setToolTip("Strain tensor elements")
		self.d11_input.setToolTip("Strain tensor elements")

		# d12 input
		self.d12_input = QLineEdit(self)
		self.d12_input.setPlaceholderText(str(self.d12))
		self.d12_input.setFixedWidth(60)
		self.d12_label = QLabel("d12", self)
		# self.d12_label.setFixedWidth(20)
		self.d12_SliderLabel = QLabel(str("%"), self) 
		# self.d12_SliderLabel.setFixedWidth(20)
		self.d12_input.returnPressed.connect(self.update_d12)
		self.d12_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.d12_btn.clicked.connect(self.update_d12)
		self.d12_btn.setAutoDefault(False)
		self.d12_label.setToolTip("Strain tensor elements")
		self.d12_input.setToolTip("Strain tensor elements")

		# d22 input
		self.d22_input = QLineEdit(self)
		self.d22_input.setPlaceholderText(str(self.d22))
		self.d22_input.setFixedWidth(60)
		self.d22_label = QLabel("d22", self)
		# self.d22_label.setFixedWidth(20)
		self.d22_SliderLabel = QLabel(str("%"), self) 
		# self.d22_SliderLabel.setFixedWidth(20)
		self.d22_input.returnPressed.connect(self.update_d22)
		self.d22_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.d22_btn.clicked.connect(self.update_d22)
		self.d22_btn.setAutoDefault(False)
		self.d22_label.setToolTip("Strain tensor elements")
		self.d22_input.setToolTip("Strain tensor elements")


		h3box = QHBoxLayout(self)
		h3box.addWidget(self.d11_label)
		h3box.addWidget(self.d11_input)
		h3box.addWidget(self.d11_SliderLabel)
		h3box.addWidget(self.d11_btn)


		h4box = QHBoxLayout(self)
		h4box.addWidget(self.d12_label)
		h4box.addWidget(self.d12_input)
		h4box.addWidget(self.d12_SliderLabel)
		h4box.addWidget(self.d12_btn)
		

		h5box = QHBoxLayout(self)
		h5box.addWidget(self.d22_label)
		h5box.addWidget(self.d22_input)
		h5box.addWidget(self.d22_SliderLabel)
		h5box.addWidget(self.d22_btn)



		# Put all the strain widgets into tab2b 
		self.tab2b.layout = QVBoxLayout(self)
		self.tab2b.layout.addWidget(self.strain2_label)
		self.tab2b.layout.addLayout(h3box)
		self.tab2b.layout.addLayout(h4box)
		self.tab2b.layout.addLayout(h5box)




		# # # # # # # # # # # # # # # # # # # # # 
		# ### Pick origin site for lattice 2 ###
		self.origin2label = QLabel("Lattice site at origin")
		# Create radiobuttons
		self.hollowsite2 = QRadioButton("Hollow")
		self.Asite2 = QRadioButton("A-site")
		self.Bsite2 = QRadioButton("B-site")

		self.hollowsite2.setChecked(True)

		# Set hover tool tips
		self.origin2label.setToolTip("Only works for hexagonal lattices.")
		self.hollowsite2.setToolTip("No atom at the origin")
		self.Asite2.setToolTip("A-atom at the origin")
		self.Bsite2.setToolTip("B-atom at the origin")
		

		# # Create QButtonGroup to be mutually exclusive  buttons. This is all u have to do for them to work correctly! 
		self.origin2_group = QButtonGroup(self)
		self.origin2_group.addButton(self.hollowsite2)
		self.origin2_group.addButton(self.Asite2)
		self.origin2_group.addButton(self.Bsite2)

		# ### IMPORTANTTTT: have to connect all buttons to the same update__button function so theyll be like mutually exclusive. so if you click one, it means the others are set to false, etc. this is defined in the updating function
		self.hollowsite2.toggled.connect(self.updateOrigin2)
		self.Asite2.toggled.connect(self.updateOrigin2)
		self.Bsite2.toggled.connect(self.updateOrigin2)

		h6box = QHBoxLayout(self)
		h6box.addWidget(self.hollowsite2)
		h6box.addWidget(self.Asite2)
		h6box.addWidget(self.Bsite2)


		# self.sublattices_label1 = QLabel("Pick strength of sublattices Z = A + B")
		# self.sublattices_label2 = QLabel("Honeycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")

		self.sublattices_label = QLabel("\n\nWeight of sublattices")
		self.sublattices_label.setToolTip("Only works if hexagonal symmetry is selected.\nHoneycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")


		# # # # # # # # # # # # # # # # # # # # 
		### Pick alpha2 weight of sublattice a ###
		# self.alpha1_param_label = QLabel("Lattice constant (nm)", self)
		# self.alpha1_param_label.setToolTip("Weight of sublattice a")
		self.alpha2_label = QLabel("alpha2:", self)
		self.alpha2_label.setToolTip("Weight of sublattice A")
		self.alpha2_input = QLineEdit(self)
		self.alpha2_input.returnPressed.connect(self.update_alpha2) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.alpha2_input.setPlaceholderText(str(self.alpha2))
		self.alpha2_input.setToolTip("Weight of sublattice A")
		self.alpha2_input.setFixedWidth(60)

		
		self.alpha2_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.alpha2_btn.clicked.connect(self.update_alpha2)
		self.alpha2_btn.setAutoDefault(False)
		self.alpha2_btn.setToolTip("Weight of sublattice A")


		h7box = QHBoxLayout(self)
		h7box.addWidget(self.alpha2_label)
		h7box.addWidget(self.alpha2_input)
		h7box.addWidget(self.alpha2_btn)



		# # # # # # # # # # # # # # # # # # # # 
		### Pick beta2 weight of sublattice b ###
		# self.beta1_param_label = QLabel("Lattice constant (nm)", self)
		# self.beta1_param_label.setToolTip("Weight of sublattice B")
		self.beta2_label = QLabel("beta2:", self)
		self.beta2_label.setToolTip("Weight of sublattice B")
		self.beta2_input = QLineEdit(self)
		self.beta2_input.returnPressed.connect(self.update_beta2) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.beta2_input.setPlaceholderText(str(self.beta2))
		self.beta2_input.setToolTip("Weight of sublattice B")
		self.beta2_input.setFixedWidth(60)

		
		self.beta2_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.beta2_btn.clicked.connect(self.update_beta2)
		self.beta2_btn.setAutoDefault(False)
		self.beta2_btn.setToolTip("Weight of sublattice B")


		h8box = QHBoxLayout(self)
		h8box.addWidget(self.beta2_label)
		h8box.addWidget(self.beta2_input)
		h8box.addWidget(self.beta2_btn)



		# Add honeycomb/dots button widgets to the tab 3
		self.tab3b.layout = QVBoxLayout(self)
		self.tab3b.layout.addWidget(self.origin2label)
		self.tab3b.layout.addLayout(h6box)
		self.tab3b.layout.addWidget(self.sublattices_label)
		self.tab3b.layout.addLayout(h7box)
		self.tab3b.layout.addLayout(h8box)

		self.tab3b.setToolTip("Only works for hexagonal lattices")


		# Set the layouts for each tab...
		self.tab1b.setLayout(self.tab1b.layout)
		self.tab1b.layout.setSpacing(1)
		self.tab2b.setLayout(self.tab2b.layout)
		self.tab2b.layout.setSpacing(1)
		self.tab3b.setLayout(self.tab3b.layout)
		self.tab3b.layout.setSpacing(1)


		vlayout = QVBoxLayout(self)
		vlayout.addWidget(self.lat2tabs)
		vlayout.setSpacing(2)

		lat2groupBox.setLayout(vlayout)
		# lat2groupBox.setContentsMargins(0,15,0,0) # Sets the left , top , right , and bottom margins to use around the layout.

		return lat2groupBox

	def initLattice3Parameters(self):
		lat3groupBox = QGroupBox("Lattice 3")


		# Add tabs to save space. 1 tab for lattice parameters (symmetry, periodicity, honeycomb/dots), the other for strain
		# https://pythonspot.com/pyqt5-tabs/
		self.lat3tabs = QTabWidget(self)

		self.tab1c = QWidget(self)
		self.tab2c = QWidget(self)
		self.tab3c = QWidget(self)
		self.lat3tabs.addTab(self.tab1c, "Parameters")
		self.lat3tabs.addTab(self.tab3c, "Sublattices")
		self.lat3tabs.addTab(self.tab2c, "Strain")


		# # # # # # # # # # # # # # # # # # # # 
		### Pick lattice 3 symmetry ###
		self.symm_label = QLabel("Symmetry:", self)
		self.symm_label.setToolTip("Choose symmetry of third lattice")
		self.hex3btn = QRadioButton("Hexagonal")
		self.sq3btn = QRadioButton("Square")
		self.hex3btn.setChecked(True)
		self.hex3btn.setToolTip("Create hexagonal lattice")
		self.sq3btn.setToolTip("Create square lattice")

		# Connect btn to update functions when clicked https://www.tutorialspoint.com/pyqt/pyqt_qpushbutton_widget.htm
		self.hex3btn.toggled.connect(self.updateLattice3)
		self.sq3btn.toggled.connect(self.updateLattice3)

		# Create QButtonGroup to be mutually exclusive w the honeycomb buttons, from https://www.programcreek.com/python/example/108083/PyQt5.QtWidgets.QButtonGroup
		self.symm3_btn_group = QButtonGroup(self)
		self.symm3_btn_group.addButton(self.hex3btn)
		self.symm3_btn_group.addButton(self.sq3btn)

		h1box = QHBoxLayout(self)
		h1box.addWidget(self.hex3btn)
		h1box.addWidget(self.sq3btn)

		# Add symmetry button widgets to tab
		self.tab1c.layout = QVBoxLayout(self)
		self.tab1c.layout.addWidget(self.symm_label)
		self.tab1c.layout.addLayout(h1box)


		# # # # # # # # # # # # # # # # # # # # 
		### Pick lattice 3 periodicity a3 ###
		self.c_param_label = QLabel("\nLattice constant (nm)", self)
		self.c_label = QLabel("c:", self)
		# self.b_label.setFixedWidth(20)
		self.c_input = QLineEdit(self)
		self.c_input.returnPressed.connect(self.update_c) # Connect this intput dialog whenever the enter/return button is pressed
		self.c_input.setPlaceholderText(str(self.c))
		self.c_input.setMinimumWidth(85)
		self.c_param_label.setToolTip("Periodicity of 3rd lattice (spacing between atoms) in nm")
		self.c_label.setToolTip("Periodicity of 3rd lattice (spacing between atoms) in nm")
		self.c_input.setToolTip("Periodicity of 3rd lattice (spacing between atoms) in nm")

		# Define label to display the value of the slider next to the slider
		self.c_nm_Label = QLabel(" nm", self) # Display the corrected value. only up to 2 decimal pts
		# self.b_nm_Label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		self.c_nm_Label.setMinimumWidth(30)
		self.c_nm_Label.setMinimumHeight(15)

		self.c_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.c_btn.clicked.connect(self.update_c)
		self.c_btn.setAutoDefault(False)
		self.c_btn.setFixedWidth(60)


		self.c_label.setMinimumHeight(15)
		self.c_input.setMinimumHeight(15)
		self.c_param_label.setMinimumHeight(15)
		self.c_btn.setMinimumHeight(15)


		h2box = QHBoxLayout(self)
		h2box.addWidget(self.c_label)
		h2box.addWidget(self.c_input)
		h2box.addWidget(self.c_nm_Label)
		h2box.addWidget(self.c_btn)

		# Add lattice constant a widgets to the tab 1
		self.tab1c.layout.addWidget(self.c_param_label)
		self.tab1c.layout.addLayout(h2box)



		


		# # # # # # # # # # # # # # # # # # # # 
		### Twist angle ###
		self.thetatw2_title = QLabel('Twist angle', self)
		self.thetatw2_input = QLineEdit(self)
		self.thetatw2_input.returnPressed.connect(self.updateThetaTw2) # Connect this intput dialog whenever the enter/return button is pressed
		self.thetatw2_input.setPlaceholderText(str(self.theta_tw2))
		self.thetatw2_input.setToolTip("Relative twist angle between 2nd & 3rd lattice")
		self.thetatw2_input.setFixedWidth(60)

		# Define label to display the value of the slider next to the slider
		self.thetatwLabel = QLabel(' deg', self)
		self.thetatwLabel.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		# self.thetatwLabel.setMinimumWidth(30)
		self.thetatwLabel.setFixedWidth(30)

		self.theta_tw2_btn = QPushButton("Twist", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.theta_tw2_btn.clicked.connect(self.updateThetaTw2)
		self.theta_tw2_btn.setAutoDefault(False)

		self.thetatw2_input.setMinimumHeight(15)
		self.thetatwLabel.setMinimumHeight(15)
		self.theta_tw2_btn.setMinimumHeight(15)

		h7box = QHBoxLayout(self)
		h7box.addWidget(self.thetatw2_input)
		h7box.addWidget(self.thetatwLabel)
		h7box.addWidget(self.theta_tw2_btn)


		# Add twist angle widgets to the tab 1
		self.tab1c.layout.addWidget(self.thetatw2_title)
		self.tab1c.layout.addLayout(h7box)




		# # # # # # # # # # # # # # # # # # # # # 
		### Pick strain on lattice 3 ###
		self.strain3_label = QLabel("Strain (in percent, i.e. 3%)")

		# f11 input
		self.f11_input = QLineEdit(self)
		self.f11_input.setPlaceholderText(str(self.f11))
		self.f11_input.setFixedWidth(60) # Shorten the width of the text box
		self.f11_label = QLabel("f11", self)
		# self.d11_label.setFixedWidth(20)
		self.f11_SliderLabel = QLabel(str("%"), self) 
		# self.d11_SliderLabel.setFixedWidth(20)
		self.f11_input.returnPressed.connect(self.update_f11)
		self.f11_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.f11_btn.clicked.connect(self.update_f11)
		self.f11_btn.setAutoDefault(False)
		self.f11_label.setToolTip("Strain tensor elements")
		self.f11_input.setToolTip("Strain tensor elements")

		# f12 input
		self.f12_input = QLineEdit(self)
		self.f12_input.setPlaceholderText(str(self.f12))
		self.f12_input.setFixedWidth(60)
		self.f12_label = QLabel("f12", self)
		# self.d12_label.setFixedWidth(20)
		self.f12_SliderLabel = QLabel(str("%"), self) 
		# self.d12_SliderLabel.setFixedWidth(20)
		self.f12_input.returnPressed.connect(self.update_f12)
		self.f12_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.f12_btn.clicked.connect(self.update_f12)
		self.f12_btn.setAutoDefault(False)
		self.f12_label.setToolTip("Strain tensor elements")
		self.f12_input.setToolTip("Strain tensor elements")

		# f22 input
		self.f22_input = QLineEdit(self)
		self.f22_input.setPlaceholderText(str(self.f22))
		self.f22_input.setFixedWidth(60)
		self.f22_label = QLabel("f22", self)
		# self.d22_label.setFixedWidth(20)
		self.f22_SliderLabel = QLabel(str("%"), self) 
		# self.d22_SliderLabel.setFixedWidth(20)
		self.f22_input.returnPressed.connect(self.update_f22)
		self.f22_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.f22_btn.clicked.connect(self.update_f22)
		self.f22_btn.setAutoDefault(False)
		self.f22_label.setToolTip("Strain tensor elements")
		self.f22_input.setToolTip("Strain tensor elements")


		h3box = QHBoxLayout(self)
		h3box.addWidget(self.f11_label)
		h3box.addWidget(self.f11_input)
		h3box.addWidget(self.f11_SliderLabel)
		h3box.addWidget(self.f11_btn)


		h4box = QHBoxLayout(self)
		h4box.addWidget(self.f12_label)
		h4box.addWidget(self.f12_input)
		h4box.addWidget(self.f12_SliderLabel)
		h4box.addWidget(self.f12_btn)
		

		h5box = QHBoxLayout(self)
		h5box.addWidget(self.f22_label)
		h5box.addWidget(self.f22_input)
		h5box.addWidget(self.f22_SliderLabel)
		h5box.addWidget(self.f22_btn)


		# Put all the strain widgets into tab2b 
		self.tab2c.layout = QVBoxLayout(self)
		self.tab2c.layout.addWidget(self.strain3_label)
		self.tab2c.layout.addLayout(h3box)
		self.tab2c.layout.addLayout(h4box)
		self.tab2c.layout.addLayout(h5box)




		# # # # # # # # # # # # # # # # # # # # # 
		# ### Pick origin site for lattice 3 ###
		self.origin3label = QLabel("Lattice site at origin")
		# Create radiobuttons
		self.hollowsite3 = QRadioButton("Hollow")
		self.Asite3 = QRadioButton("A-site")
		self.Bsite3 = QRadioButton("B-site")

		self.hollowsite3.setChecked(True)

		# Set hover tool tips
		self.origin3label.setToolTip("Only works for hexagonal lattices.")
		self.hollowsite3.setToolTip("No atom at the origin")
		self.Asite3.setToolTip("A-atom at the origin")
		self.Bsite3.setToolTip("B-atom at the origin")
		

		# # Create QButtonGroup to be mutually exclusive  buttons. This is all u have to do for them to work correctly! 
		self.origin3_group = QButtonGroup(self)
		self.origin3_group.addButton(self.hollowsite3)
		self.origin3_group.addButton(self.Asite3)
		self.origin3_group.addButton(self.Bsite3)

		# ### IMPORTANTTTT: have to connect all buttons to the same update__button function so theyll be like mutually exclusive. so if you click one, it means the others are set to false, etc. this is defined in the updating function
		self.hollowsite3.toggled.connect(self.updateOrigin3)
		self.Asite3.toggled.connect(self.updateOrigin3)
		self.Bsite3.toggled.connect(self.updateOrigin3)

		h6box = QHBoxLayout(self)
		h6box.addWidget(self.hollowsite3)
		h6box.addWidget(self.Asite3)
		h6box.addWidget(self.Bsite3)


		# self.sublattices_label1 = QLabel("Pick strength of sublattices Z = A + B")
		# self.sublattices_label2 = QLabel("Honeycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")

		self.sublattices_label = QLabel("\n\nWeight of sublattices")
		self.sublattices_label.setToolTip("Only works if hexagonal symmetry is selected.\nHoneycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")


	

		self.sublattices_label1 = QLabel("Pick strength of sublattices (Z = A + B)")
		self.sublattices_label2 = QLabel("Honeycomb lattice: alpha = beta \nTriangular lattice: alpha = 1, beta = 0")
		# # # # # # # # # # # # # # # # # # # # 
		### Pick alpha3 weight of sublattice a ###
		# self.alpha1_param_label = QLabel("Lattice constant (nm)", self)
		# self.alpha1_param_label.setToolTip("Weight of sublattice a")
		self.alpha3_label = QLabel("alpha3:", self)
		self.alpha3_label.setToolTip("Weight of sublattice A")
		self.alpha3_input = QLineEdit(self)
		self.alpha3_input.returnPressed.connect(self.update_alpha3) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.alpha3_input.setPlaceholderText(str(self.alpha3))
		self.alpha3_input.setToolTip("Weight of sublattice A")
		self.alpha3_input.setFixedWidth(60)

		
		self.alpha3_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.alpha3_btn.clicked.connect(self.update_alpha3)
		self.alpha3_btn.setAutoDefault(False)
		self.alpha3_btn.setToolTip("Weight of sublattice A")


		h7box = QHBoxLayout(self)
		h7box.addWidget(self.alpha3_label)
		h7box.addWidget(self.alpha3_input)
		h7box.addWidget(self.alpha3_btn)



		# # # # # # # # # # # # # # # # # # # # 
		### Pick beta1 weight of sublattice b ###
		# self.beta1_param_label = QLabel("Lattice constant (nm)", self)
		# self.beta1_param_label.setToolTip("Weight of sublattice B")
		self.beta3_label = QLabel("beta3:", self)
		self.beta3_label.setToolTip("Weight of sublattice B")
		self.beta3_input = QLineEdit(self)
		self.beta3_input.returnPressed.connect(self.update_beta3) # Connect this intput dialog whenever the enter/return/tab button is pressed or you click away from the widget box
		self.beta3_input.setPlaceholderText(str(self.beta3))
		self.beta3_input.setToolTip("Weight of sublattice B")
		self.beta3_input.setFixedWidth(60)

		
		self.beta3_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update! connect to the same update function!
		self.beta3_btn.clicked.connect(self.update_beta3)
		self.beta3_btn.setAutoDefault(False)
		self.beta3_btn.setToolTip("Weight of sublattice B")


		h8box = QHBoxLayout(self)
		h8box.addWidget(self.beta3_label)
		h8box.addWidget(self.beta3_input)
		h8box.addWidget(self.beta3_btn)






		

		# Add honeycomb/dots button widgets to the tab 3
		self.tab3c.layout = QVBoxLayout(self)
		self.tab3c.layout.addWidget(self.origin3label)
		self.tab3c.layout.addLayout(h6box)
		self.tab3c.layout.addWidget(self.sublattices_label)
		self.tab3c.layout.addLayout(h7box)
		self.tab3c.layout.addLayout(h8box)

		self.tab3c.setToolTip("Only works for hexagonal lattices")




		# Set the layouts for each tab...
		self.tab1c.setLayout(self.tab1c.layout)
		self.tab1c.layout.setSpacing(1)
		self.tab2c.setLayout(self.tab2c.layout)
		self.tab2c.layout.setSpacing(1)
		self.tab3c.setLayout(self.tab3c.layout)
		self.tab3c.layout.setSpacing(1)


		# Add tabs to a QVBoxLayout to be able to set the whole groupbox layout to the tabs
		vlayout = QVBoxLayout(self)
		vlayout.addWidget(self.lat3tabs)
		vlayout.setSpacing(2)

		lat3groupBox.setLayout(vlayout)
		# lat3groupBox.setContentsMargins(0,15,0,0) # Sets the left , top , right , and bottom margins to use around the layout.

		return lat3groupBox

	def updateLattice1(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == "Hexagonal":
				self.lattice1 = "Hexagonal"
			elif radio_btn.text() == "Square":
				self.lattice1 = "Square"
		 # Will update lattice1 to "Hexagonal" or "Square" depending on which is clicked
		self.plotAtoms()

	def updateLattice2(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == "Hexagonal":
				self.lattice2 = "Hexagonal"
			elif radio_btn.text() == "Square":
				self.lattice2 = "Square"
		 # Will update lattice2 to "Hexagonal" or "Square" depending on which is clicked
		self.plotAtoms()

	def updateLattice3(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == "Hexagonal":
				self.lattice3 = "Hexagonal"
			elif radio_btn.text() == "Square":
				self.lattice3 = "Square"
		 # Will update lattice3 to "Hexagonal" or "Square" depending on which is clicked
		self.plotAtoms()

	def update_a(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.a = eval(self.a_input.text()) #float(self.a_input.text())
			self.a_nm_Label.setText('nm') 
			self.a_input.setPlaceholderText(str(self.a))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.a_error = QMessageBox()
			self.a_error.setWindowTitle("Error")
			self.a_error.setText("Type in a number or numerical expression")
			self.a_error.setInformativeText("Your input for a is: " +  str(self.a_input.text()))
			self.a_error.setIcon(QMessageBox.Warning)
			self.a_error.setStandardButtons(QMessageBox.Retry)
			x = self.a_error.exec()
					
	def update_b(self):
		try: 
			self.b = eval(self.b_input.text())
			self.b_nm_Label.setText('nm')
			self.b_input.setPlaceholderText(str(self.b))
			self.plotAtoms() 
		except:
			self.b_error = QMessageBox()
			self.b_error.setWindowTitle("Error")
			self.b_error.setText("Type in a number or numerical expression")
			self.b_error.setInformativeText("Your input for b is: " +  str(self.b_input.text()))
			self.b_error.setIcon(QMessageBox.Warning)
			self.b_error.setStandardButtons(QMessageBox.Retry)
			x = self.b_error.exec()

	def update_c(self):
		try: 
			self.c = eval(self.c_input.text())
			self.c_nm_Label.setText('nm')
			self.c_input.setPlaceholderText(str(self.c))
			self.plotAtoms() 
		except:
			self.c_error = QMessageBox()
			self.c_error.setWindowTitle("Error")
			self.c_error.setText("Type in a number or numerical expression")
			self.c_error.setInformativeText("Your input for c is: " +  str(self.c_input.text()))
			self.c_error.setIcon(QMessageBox.Warning)
			self.c_error.setStandardButtons(QMessageBox.Retry)
			x = self.c_error.exec()

	def updateThetaTw(self):
		try:
			self.theta_tw = eval(self.thetatw_input.text())
			self.thetatwLabel.adjustSize()
			self.thetatw_input.setPlaceholderText(str(self.theta_tw))
			self.plotAtoms() 
		except:
			self.theta_tw_error = QMessageBox()
			self.theta_tw_error.setWindowTitle("Error")
			self.theta_tw_error.setText("Type in a number or numerical expression")
			self.theta_tw_error.setInformativeText("Your input for theta_twist is: " +  str(self.thetatw_input.text()))
			self.theta_tw_error.setIcon(QMessageBox.Warning)
			self.theta_tw_error.setStandardButtons(QMessageBox.Retry)
			x = self.theta_tw_error.exec()

	def updateThetaTw2(self):
		try:
			self.theta_tw2 = eval(self.thetatw2_input.text())
			self.thetatwLabel.adjustSize()
			self.thetatw2_input.setPlaceholderText(str(self.theta_tw2))
			self.plotAtoms() 
		except:
			self.theta_tw2_error = QMessageBox()
			self.theta_tw2_error.setWindowTitle("Error")
			self.theta_tw2_error.setText("Type in a number or numerical expression")
			self.theta_tw2_error.setInformativeText("Your input for theta_twist2 is: " +  str(self.thetatw2_input.text()))
			self.theta_tw2_error.setIcon(QMessageBox.Warning)
			self.theta_tw2_error.setStandardButtons(QMessageBox.Retry)
			x = self.theta_tw2_error.exec()

	def update_e11(self):
		try: 
			self.e11 = eval(self.e11_input.text())
			self.e11 /= 100 # Divide by 100 to scale the strain down to percent
			self.e11_input.setPlaceholderText(str(self.e11))
			self.plotAtoms() 
		except:
			self.e11_error = QMessageBox()
			self.e11_error.setWindowTitle("Error")
			self.e11_error.setText("Type in a number or numerical expression")
			self.e11_error.setInformativeText("Your input for e11 is: " +  str(self.e11_input.text()))
			self.e11_error.setIcon(QMessageBox.Warning)
			self.e11_error.setStandardButtons(QMessageBox.Retry)
			x = self.e11_error.exec()

	def update_e12(self):
		try:
			self.e12 = eval(self.e12_input.text())
			self.e12 /=100 # Divide by 100 to scale the strain down to percent
			self.e12_input.setPlaceholderText(str(self.e12))
			self.plotAtoms()
		except:
			self.e12_error = QMessageBox()
			self.e12_error.setWindowTitle("Error")
			self.e12_error.setText("Type in a number or numerical expression")
			self.e12_error.setInformativeText("Your input for e12 is: " +  str(self.e12_input.text()))
			self.e12_error.setIcon(QMessageBox.Warning)
			self.e12_error.setStandardButtons(QMessageBox.Retry)
			x = self.e12_error.exec()

	def update_e22(self):
		try: 
			self.e22 = eval(self.e22_input.text())
			self.e22_input.setPlaceholderText(str(self.e22))
			self.e22 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms()
		except:
			self.e22_error = QMessageBox()
			self.e22_error.setWindowTitle("Error")
			self.e22_error.setText("Type in a number or numerical expression")
			self.e22_error.setInformativeText("Your input for e22 is: " +  str(self.e22_input.text()))
			self.e22_error.setIcon(QMessageBox.Warning)
			self.e22_error.setStandardButtons(QMessageBox.Retry)
			x = self.e22_error.exec()

	def update_d11(self):
		try:
			self.d11 = eval(self.d11_input.text())
			self.d11_input.setPlaceholderText(str(self.d11))
			self.d11 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms() 
		except:
			self.d11_error = QMessageBox()
			self.d11_error.setWindowTitle("Error")
			self.d11_error.setText("Type in a number or numerical expression")
			self.d11_error.setInformativeText("Your input for d11 is: " +  str(self.d11_input.text()))
			self.d11_error.setIcon(QMessageBox.Warning)
			self.d11_error.setStandardButtons(QMessageBox.Retry)
			x = self.d11_error.exec()

	def update_d12(self):
		try:
			self.d12 = eval(self.d12_input.text())
			self.d12_input.setPlaceholderText(str(self.d12))
			self.d12 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms() 
		except:
			self.d12_error = QMessageBox()
			self.d12_error.setWindowTitle("Error")
			self.d12_error.setText("Type in a number or numerical expression")
			self.d12_error.setInformativeText("Your input for d12 is: " +  str(self.d12_input.text()))
			self.d12_error.setIcon(QMessageBox.Warning)
			self.d12_error.setStandardButtons(QMessageBox.Retry)
			x = self.d12_error.exec()

	def update_d22(self):
		try:
			self.d22 = eval(self.d22_input.text())
			self.d22_input.setPlaceholderText(str(self.d22))
			self.d22 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms() 
		except:
			self.d22_error = QMessageBox()
			self.d22_error.setWindowTitle("Error")
			self.d22_error.setText("Type in a number or numerical expression")
			self.d22_error.setInformativeText("Your input for d22 is: " +  str(self.d22_input.text()))
			self.d22_error.setIcon(QMessageBox.Warning)
			self.d22_error.setStandardButtons(QMessageBox.Retry)
			x = self.d22_error.exec()
	
	def update_f11(self):
		try:
			self.f11 = eval(self.f11_input.text())
			self.f11_input.setPlaceholderText(str(self.f11))
			self.f11 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms() 
		except:
			self.f11_error = QMessageBox()
			self.f11_error.setWindowTitle("Error")
			self.f11_error.setText("Type in a number or numerical expression")
			self.f11_error.setInformativeText("Your input for d11 is: " +  str(self.f11_input.text()))
			self.f11_error.setIcon(QMessageBox.Warning)
			self.f11_error.setStandardButtons(QMessageBox.Retry)
			x = self.f11_error.exec()

	def update_f12(self):
		try:
			self.f12 = eval(self.f12_input.text())
			self.f12_input.setPlaceholderText(str(self.f12))
			self.f12 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms() 
		except:
			self.f12_error = QMessageBox()
			self.f12_error.setWindowTitle("Error")
			self.f12_error.setText("Type in a number or numerical expression")
			self.f12_error.setInformativeText("Your input for f12 is: " +  str(self.f12_input.text()))
			self.f12_error.setIcon(QMessageBox.Warning)
			self.f12_error.setStandardButtons(QMessageBox.Retry)
			x = self.f12_error.exec()

	def update_f22(self):
		try:
			self.f22 = eval(self.f22_input.text())
			self.f22_input.setPlaceholderText(str(self.f22))
			self.f22 /= 100 # Divide by 100 to scale the strain down to percent
			self.plotAtoms() 
		except:
			self.f22_error = QMessageBox()
			self.f22_error.setWindowTitle("Error")
			self.f22_error.setText("Type in a number or numerical expression")
			self.f22_error.setInformativeText("Your input for f22 is: " +  str(self.f22_input.text()))
			self.f22_error.setIcon(QMessageBox.Warning)
			self.f22_error.setStandardButtons(QMessageBox.Retry)
			x = self.f22_error.exec()

	def updateOrigin1(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == 'Hollow':
				self.origin1 = "Hollow"
			elif radio_btn.text() == 'A-site':
				self.origin1 = 'A-site'
			elif radio_btn.text() == "B-site":
				self.origin1 = "B-site"
		self.plotAtoms() 

	def updateOrigin2(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == 'Hollow':
				self.origin2 = "Hollow"
			elif radio_btn.text() == 'A-site':
				self.origin2 = 'A-site'
			elif radio_btn.text() == "B-site":
				self.origin2 = "B-site"

		self.plotAtoms() 

	def updateOrigin3(self):
		radio_btn = self.sender()
		if radio_btn.isChecked():
			if radio_btn.text() == 'Hollow':
				self.origin3 = "Hollow"
			elif radio_btn.text() == 'A-site':
				self.origin3 = 'A-site'
			elif radio_btn.text() == "B-site":
				self.origin3 = "B-site"
		self.plotAtoms() 

	def update_alpha1(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.alpha1 = eval(self.alpha1_input.text()) 
			self.alpha1_input.setPlaceholderText(str(self.alpha1))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.alpha1_error = QMessageBox()
			self.alpha1_error.setWindowTitle("Error")
			self.alpha1_error.setText("Type in a number or numerical expression")
			self.alpha1_error.setInformativeText("Your input for alpha1 is: " +  str(self.alpha1_input.text()))
			self.alpha1_error.setIcon(QMessageBox.Warning)
			self.alpha1_error.setStandardButtons(QMessageBox.Retry)
			x = self.alpha1_error.exec()

	def update_beta1(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.beta1 = eval(self.beta1_input.text()) 
			self.beta1_input.setPlaceholderText(str(self.beta1))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.beta1_error = QMessageBox()
			self.beta1_error.setWindowTitle("Error")
			self.beta1_error.setText("Type in a number or numerical expression")
			self.beta1_error.setInformativeText("Your input for beta1 is: " +  str(self.beta1_input.text()))
			self.beta1_error.setIcon(QMessageBox.Warning)
			self.beta1_error.setStandardButtons(QMessageBox.Retry)
			x = self.beta1_error.exec()

	def update_alpha2(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.alpha2 = eval(self.alpha2_input.text()) 
			self.alpha2_input.setPlaceholderText(str(self.alpha2))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.alpha2_error = QMessageBox()
			self.alpha2_error.setWindowTitle("Error")
			self.alpha2_error.setText("Type in a number or numerical expression")
			self.alpha2_error.setInformativeText("Your input for alpha2 is: " +  str(self.alpha2_input.text()))
			self.alpha2_error.setIcon(QMessageBox.Warning)
			self.alpha2_error.setStandardButtons(QMessageBox.Retry)
			x = self.alpha2_error.exec()

	def update_beta2(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.beta2 = eval(self.beta2_input.text()) 
			self.beta2_input.setPlaceholderText(str(self.beta2))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.beta2_error = QMessageBox()
			self.beta2_error.setWindowTitle("Error")
			self.beta2_error.setText("Type in a number or numerical expression")
			self.beta2_error.setInformativeText("Your input for beta2 is: " +  str(self.beta2_input.text()))
			self.beta2_error.setIcon(QMessageBox.Warning)
			self.beta2_error.setStandardButtons(QMessageBox.Retry)
			x = self.beta2_error.exec()

	def update_alpha3(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.alpha3 = eval(self.alpha3_input.text()) 
			self.alpha3_input.setPlaceholderText(str(self.alpha3))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.alpha3_error = QMessageBox()
			self.alpha3_error.setWindowTitle("Error")
			self.alpha3_error.setText("Type in a number or numerical expression")
			self.alpha3_error.setInformativeText("Your input for alpha3 is: " +  str(self.alpha3_input.text()))
			self.alpha3_error.setIcon(QMessageBox.Warning)
			self.alpha3_error.setStandardButtons(QMessageBox.Retry)
			x = self.alpha3_error.exec()

	def update_beta3(self):
		try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
				# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			self.beta3 = eval(self.beta3_input.text()) 
			self.beta3_input.setPlaceholderText(str(self.beta3))
			self.plotAtoms() 
		except:
			# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# https://www.w3schools.com/python/python_try_except.asp
			# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			self.beta3_error = QMessageBox()
			self.beta3_error.setWindowTitle("Error")
			self.beta3_error.setText("Type in a number or numerical expression")
			self.beta3_error.setInformativeText("Your input for beta3 is: " +  str(self.beta3_input.text()))
			self.beta3_error.setIcon(QMessageBox.Warning)
			self.beta3_error.setStandardButtons(QMessageBox.Retry)
			x = self.beta3_error.exec()

	# CREATING THE MATPLOTLIB FIGURE
	def initMatplotlibFig(self):
		# Create a figure instance to plot on
		self.figure = plt.figure(figsize=(10,10))

		# Create a grid to plot mutiple plots on
		self.grid = GridSpec(nrows=1,ncols=2)

		groupbox = QGroupBox()

		# creating a Vertical Box layout for storing widgets
		vlayout = QVBoxLayout()

		# this is the Canvas Widget that displays the 'figure'
		# it takes the 'figure' instance as a parameter to __init__
		self.canvas = FigureCanvas(self.figure)

		# This connects a button press event when u click on a plot, to the class method self.onClickEvent
		# When this is activated it performs whatever actions are specified in self.onClickEvent functin
		# self.cid = self.figure.canvas.mpl_connect('button_press_event', self.onClickEvent)

		# this is the Navigation widget for matplotlib figs - includes zoom/pan/save/other cool feats
		# it takes the Canvas widget and a parent
		self.toolbar = NavigationToolbar(self.canvas, self)
		# self.toolbar.setStyleSheet("color: red;background: black")
		# self.toolbar.setStyleSheet("background-color:magenta;")
		# groupbox.setStyleSheet("background-color:magenta;")
		# self.canvas.setStyleSheet("background-color:magenta;")
		# adding tool bar to the layout
		vlayout.addWidget(self.toolbar)

		# adding canvas to the layout
		vlayout.addWidget(self.canvas)

		groupbox.setLayout(vlayout)

		plt.tight_layout()

		self.plotAtoms()

		return groupbox

	# COLORBAR FUNCTION TO FIX COLORBar issues. stolen from https://joseph-long.com/writing/colorbars/
	# not using this in the gui anymore bc it caused issues. leaving it here for future reference only
	def make_colorbar(self, mappable):
		from mpl_toolkits.axes_grid1 import make_axes_locatable
		import matplotlib.pyplot as plt
		last_axes = plt.gca()
		ax = mappable.axes
		fig = ax.figure
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = fig.colorbar(mappable, cax=cax)
		plt.sca(last_axes)
		return cbar

	def initColormapDropdown(self):
		groupBox = QGroupBox("Colormap")

		# Add a dropdown list to choose colormap to use for plotting the REAL SPACE image
		self.dropdownColormap_RS = QComboBox(self)
		self.dropdownColormap_RS.addItems(self.colormapList)
		self.dropdownColormap_RS.setCurrentText(self.colormap_RS) 

		# This is/updates the label/text to the chosen text from the dropdown/combo list
		self.colormapLabel_RS= QLabel("Real space image: ", self) # idk what this actually is for

		# This sets the current value of colormap instance variable (the one thats actually used in the plotting functions) to the default / current value
		self.colormap_RS = self.dropdownColormap_RS.currentText()

		# To update which text appears in the coombo box to the chosen text from the dropdown/combo list
		self.dropdownColormap_RS.activated[str].connect(self.updateColormap)



		# Add a dropdown list to choose colormap to use for plotting the REAL SPACE image
		self.dropdownColormap_FFT = QComboBox(self)
		self.dropdownColormap_FFT.addItems(self.colormapList)
		self.dropdownColormap_FFT.setCurrentText(self.colormap_FFT)

		# This is/updates the label/text to the chosen text from the dropdown/combo list
		self.colormapLabel_FFT = QLabel("FFT:", self) # idk what this actually is for

		# This sets the current value of colormap instance variable (the one thats actually used in the plotting functions) to the default / current value
		self.colormap_FFT = self.dropdownColormap_FFT.currentText()

		# To update which text appears in the coombo box to the chosen text from the dropdown/combo list
		self.dropdownColormap_FFT.activated[str].connect(self.updateColormap)


		vlayout = QVBoxLayout()


		# Place the combobox
		hlayout = QHBoxLayout()
		hlayout.addWidget(self.colormapLabel_RS)
		hlayout.addWidget(self.dropdownColormap_RS)


		hlayout2 = QHBoxLayout()
		hlayout2.addWidget(self.colormapLabel_FFT)
		hlayout2.addWidget(self.dropdownColormap_FFT)

		vlayout.addLayout(hlayout)
		vlayout.addLayout(hlayout2)

		groupBox.setLayout(vlayout) # Add the widgets and value text to the groupbox

		return groupBox

	def updateColormap(self):
		# self.colormap_RS.setText(text)
		# This updates the value of colormap whenever its changed and so when you click the buttoon it updates the plot to the selected value from dropdown list
		self.colormap_RS = self.dropdownColormap_RS.currentText()
		self.colormap_FFT = self.dropdownColormap_FFT.currentText()
		self.plotAtoms()

	# THIS IS THE ACTUAL PLOTTING FUNCTION
	def plotAtoms(self):
		# clearing old figure
		self.figure.clear()

		# if moire btn is clicked yes, plot a moire lattice
		if self.moireBtn != 'Single':
			self.Z, self.fftZ = moirelattice(self.pix, self.L, self.a, self.b, self.c, self.moireBtn, self.lattice1, self.lattice2, self.lattice3, self.theta_im, self.theta_tw, self.theta_tw2, self.e11, self.e12, self.e22, self.d11, self.d12, self.d22, self.f11, self.f12, self.f22, self.alpha1, self.beta1, self.alpha2, self.beta2, self.alpha3, self.beta3, self.origin1, self.origin2, self.origin3, self.filter_bool, self.sigma)

		# if moire btn is clicked no, only plot a single lattice using the lattice1 parameters. all lattice2 inputs are ignored
		else:
			if self.lattice1 == 'Hexagonal':
				self.Z, self.fftZ = hexatoms(self.pix, self.L, self.a, self.theta_im, self.e11, self.e12, self.e22, self.alpha1, self.beta1, self.origin1)
			elif self.lattice1 == 'Square':
				self.Z, self.fftZ = squareatoms(self.pix, self.L, self.a, self.theta_im, self.e11, self.e12, self.e22)

			# Normalize the FFTs to be between 0-1 (bc hexatoms only normalizes Z, moirelattice is what normalizes fftZ, but if you chose single lattice, moirelattice code isnt run. so need to normalize the FFT here)
			self.fftZ = (self.fftZ - np.min(np.min((self.fftZ))))/(np.max(np.max(self.fftZ)) - np.min(np.min(self.fftZ)))
	 


		# Create an axis for plotting - the matplotlib gridspec figure was defined in self.__init__ in the GUI app
		ax00 = self.figure.add_subplot(self.grid[0,0])
		ax01 = self.figure.add_subplot(self.grid[0,1])



		# Fix to correct for issue with Matplotlib version 3.5.1 creating cartoonishly large fonts
		import matplotlib as mpl
		mpl_ver = matplotlib.__version__.split('.')

		if mpl_ver[0] == 3 and mpl_ver[1] == 5:
			# Set new font sizes and spacings
			plt.rcParams.update({'font.size': 5})



		### FYI: I do np.flipud(image) and then ax.invert_axes to set the y-axis to be normal like from 0--->N strting at the origin instead of it being reversed, from N ---> 0
		# Plot original image
		fig1 = ax00.imshow(np.flipud(self.Z), cmap = self.colormap_RS, extent=[0,self.L,0,self.L])
		ax00.set_xticks([])
		ax00.set_yticks([])
		ax00.set_title('Real space image')
		ax00.grid(False)
		plt.colorbar(fig1, ax=ax00, fraction=0.046, pad=0.04) # fixed colorbar issues, from: https://stackoverflow.com/questions/16702479/matplotlib-colorbar-placement-and-size
		plt.tight_layout()

		# Plot a circular with the radius of the half-width at half-max of a 2D gaussian of width w, HWHM = sqrt(2*log(2))*w
		if self.filter_bool == True and self.sigma != 0:
			self.sigma_real = self.sigma*self.L/(self.pix-1) # Convert pixels to real-space units
			circ = plt.Circle((self.L/8,self.L/8),self.sigma_real*np.sqrt(2*np.log(2)),fill=True,color='white',alpha=0.5)
			ax00.add_artist(circ)

		from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
		import matplotlib.font_manager as fm
		# fontprops = fm.FontProperties(size=12)
		scalebar_Z = AnchoredSizeBar(ax00.transData, np.ceil(self.L/10 *2)/2, str(np.ceil(self.L/10 *2)/2) + ' nm', 
					     'lower right', pad=0.5, sep=5, borderpad=0.5, frameon=True, color='black', label_top=True)

		ax00.add_artist(scalebar_Z)

		# Plot FFT of original image
		
		# Calculate min/max extent of k-space given N pixels and L distance for a sampling rate L/N
		freq = np.fft.fftfreq(self.pix, self.L/self.pix)
		extL = 2*np.pi*np.min(freq) 
		extR = 2*np.pi*np.max(freq)
		
		# # Below is an explicit calculation to see how Numpy sorts the frequencies for even/odd N pixels
		# if self.pix%2==0:
		# 	extL = -self.pix*np.pi/self.L
		# 	extR = (self.pix-2)*np.pi/self.L
		# else:
		# 	extL = -(self.pix-1)*np.pi/self.L
		# 	extR = (self.pix-1)*np.pi/self.L
		
		fig2 = ax01.imshow(np.flipud(self.fftZ), cmap = self.colormap_FFT, extent=[extL, extR, extL, extR])
		ax01.set_xlabel('$k_x$ (1/nm)')
		ax01.set_ylabel('$k_y$ (1/nm)', labelpad= -20)
		ax01.set_title('FFT')
		ax01.grid(False)

		plt.colorbar(fig2, ax=ax01, fraction=0.046, pad=0.04)

		plt.tight_layout()

		if self.filter_bool == True and self.sigma != 0:
			sigma_k = 1/self.sigma_real # The width of gaussian in k-space
			circ_k = plt.Circle((0,0),sigma_k*np.sqrt(2*np.log(2)),fill=False,color='red')
			ax01.add_artist(circ_k)

		# scalebar_FFT = AnchoredSizeBar(ax01.transData, np.ceil(np.pi/(4*self.a)), str(np.ceil(np.pi/self.a)) + r' $\mathrm{nm}^{-1}$',
		# 			       'lower right', pad=0.5, sep=5, borderpad=0.5, frameon=True, color='black',label_top=True)

		# ax01.add_artist(scalebar_FFT)
		self.canvas.draw() 

	def initSaveButton(self):
		groupBox = QGroupBox("Save files")
		groupBox.setToolTip("Creates a new folder and saves the real space image array as a .txt file,\na new .txt file with the input parameters\nand .png files of each image")
		
		# # Add Qline Edit to enter file name
		# self.saveFileName = QLineEdit(self)
		# self.saveFileName.returnPressed.connect(self.updateSaveButton) # Connect this intput dialog whenever the enter/return button is pressed
		# self.saveFileName.setPlaceholderText('Save as...')

		self.saveLabel = QLabel("Input file name only (no extension)")
		self.saveLabel.setWordWrap(True)
		self.saveLabel.setMinimumHeight(20)
		# Add QPushButton to open file directory whenever you want to save
		self.save_btn = QPushButton("Click to save", self)
		self.save_btn.clicked.connect(self.updateSaveButton)
		self.save_btn.setAutoDefault(False) # This is so its not default on, you can press enter on other widgets and this window wont pop up
		
		vlayout = QVBoxLayout()
		vlayout.addWidget(self.saveLabel)
		vlayout.addWidget(self.save_btn)
		groupBox.setLayout(vlayout)

		return groupBox

	def updateSaveButton(self):
		# The line below opens the file dirctory so you can choose where to save the files...
		# from https://www.tutorialexample.com/an-introduce-to-pyqt-qfiledialog-get-directory-path-with-examples-pyqt-tutorial/
		# & https://www.tutorialexample.com/an-introduce-to-pyqt-button-bind-click-event-with-examples-pyqt-tutorial/ 
		# & https://pythonspot.com/pyqt5-file-dialog/
		# self.dir_path = QFileDialog.getExistingDirectory(self, "Choose Directory", os.getcwd()) # This just opens the directory files. this seems like its more for opening/selecting a certain file rather than saving it 
		self.saveFileName = QFileDialog.getSaveFileName(self, "Save file as...",os.getcwd())  # This returns a 2 item tuple, and the first item is the path directory w/the file name so save the files here
		# Use try/except because the gui would crash if i opened the file dialog to save a file, then pressed cancel
		try:
			self.filePath=self.saveFileName[0] # This includes the full path... has form: ~/users/Desktop/folder/filename
			self.fileName = self.filePath.split('/')[-1] # Split the string by the slant / character. Splitting it creates a list, the [-1] index retrieves the LAST entry of that list, which is going to be the main file name of all the files that will be saved
															

			# Check whether the specified path exists or not. from https://appdividend.com/2021/07/03/how-to-create-directory-if-not-exist-in-python/
			isExist = os.path.exists(self.filePath)
			if not isExist:
				# Create a new directory because it does not exist 
				os.makedirs(self.filePath) # of the form ~/users/Desktop/filename (file name has no extension!)

		

			# Save lattice & FFT figure as png's in the created folder:
			# plt.savefig(self.fileName+'.png')
			mplimg.imsave(self.filePath + '/' + self.fileName+'.png', self.Z, cmap = self.colormap_RS)
			mplimg.imsave(self.filePath + '/' + self.fileName+'_FFT.png', self.fftZ, cmap = self.colormap_FFT)
			
			# Save atoms array as txt file in the created folder
			np.savetxt(self.filePath + '/' + self.fileName +'.txt', self.Z)

			# Also save a .txt file with the input parameter values: (IN THE FUTURE maybe make this a pandas dataframe or csv file or something?)
			param_file = open(self.filePath + '/' + self.fileName +'_params.txt', "w+") # Open a new blank text file where we will write the input parameters
			param_file.write("# of layers: " + (self.moireBtn) + "\nPix: " + str(self.pix) + "\nL: " + str(self.L) +  "\nImage offset angle: " + str(self.theta_im) + "\nLow pass filter: " + str(self.filter_bool) + '\nSigma: ' + str(self.sigma) +
					'\n\n--------------------------------\nLattice 1:\n--------------------------------\n' + self.lattice1 + '\na: ' + str(self.a) + '\ne11: ' + str(self.e11) +
					'\ne12: ' + str(self.e12) + '\ne22: ' + str(self.e22) + "\nAlpha1: " + str(self.alpha1) + "\nBeta1: " + str(self.beta1) + "\nOrigin1: " + str(self.origin1) + 
					'\n\n--------------------------------\nLattice 2:\n--------------------------------\n' + self.lattice2 + '\nb: ' + str(self.b) + '\nd11: ' + str(self.d11) + 
					'\nd12: ' + str(self.d12) + '\nd22: ' + str(self.d22) + "\nAlpha2: " + str(self.alpha2) + "\nBeta2: " + str(self.beta2) + "\nOrigin2: " + str(self.origin2) + 
					'\nTwist angle (btwn lattice 1 & 2): ' + str(self.theta_tw) + '\n\n--------------------------------\nLattice 3:\n--------------------------------\n' + self.lattice3 + '\nc: ' + str(self.c) + '\nf11: ' + str(self.f11) + 
					'\nf12: ' + str(self.f12) + '\nf22: ' + str(self.f22) + "\nAlpha3: " + str(self.alpha3) + "\nBeta3: " + str(self.beta3) + "\nOrigin3: " + str(self.origin3) + 
					'\nTwist angle (btwn lattice 2 & 3): ' + str(self.theta_tw2) + 
					'\n\n--------------------------------\nEstimated time to take STM topography: ' + str(self.hrs) + 'h ' + str(self.mins) + 'min ' + str(self.sec) + ' s\nTip scanner speed: ' + str(self.vt) + ' nm/s'
					'\nEstimated time to take dI/dV map: ' + str(self.days_map) + 'days ' + str(self.hrs_map) + 'h ' + str(self.mins_map) +  'min\nTime per spectra: ' + str(self.tps) + ' s')
					
			param_file.close() # Close the text file
		

		except FileNotFoundError: # https://www.pythonfixing.com/2022/03/fixed-pyqt5-application-shutting-down.html
								  # If you press cancel when the file dialog window opens, it will return self.saveFileName = '' (an empty string), which gives error in the rest of the code above
			return 				  # ^^ so it will bring up a FileNotFound error; if this happens, just do nothing / the filedialog box will close but the app wont crash 

	def initFiltering(self):
		groupBox = QGroupBox("Low pass filtering")
		groupBox.setToolTip("Low pass filter the data")
		vlayout = QVBoxLayout(self)

		self.filter_btn = QCheckBox("Filter (pix)")
		self.filter_btn.setChecked(False)
		self.filter_btn.stateChanged.connect(self.updateSigma)
		self.filter_btn.setToolTip("Check to filter")
		# self.calc_btn_grp = QButtonGroup()
		# self.calc_btn_grp.addButton(self.calc_topo_btn, 1)

		

		# Create vt input text box
		self.sigma_input = QLineEdit(self)
		self.sigma_input.returnPressed.connect(self.updateSigma) # Connect this intput dialog whenever the enter/return button is pressed
		self.sigma_input.setPlaceholderText(str(self.sigma))
		self.sigma_input.setFixedWidth(60)
		self.sigma_input.setToolTip("Set size of gaussian mask (in real space)")

		self.sigma_btn = QPushButton("Go", self) # Create a QPushButton so users can press enter and/or click this button to update!
		self.sigma_btn.clicked.connect(self.updateSigma)
		self.sigma_btn.setAutoDefault(False)

		# Define label to display the value of the slider next to the textbox
		self.sigma_label = QLabel("Sigma:", self)

		self.sigma_real_label = QLabel("     %.2f nm"  % (self.sigma_real), self)
		# self.sigma_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		# self.sigma_label.setMinimumWidth(40)
		# self.vt_nm_label = QLabel(" nm/s", self)
		# self.vt_nm_label.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
		# self.vt_nm_label.setMinimumWidth(30)


		hbox = QHBoxLayout(self)
		hbox.addWidget(self.filter_btn)
		hbox.addWidget(self.sigma_real_label)

		vlayout = QVBoxLayout(self)
		# vlayout.addLayout(hbox)


		hbox2 = QHBoxLayout(self)
		hbox2.addWidget(self.sigma_label)
		hbox2.addWidget(self.sigma_input)
		# hbox.addWidget(self.vt_nm_label)
		hbox2.addWidget(self.sigma_btn)
		# hbox2.setSpacing(2)

		vlayout.addLayout(hbox)
		vlayout.addLayout(hbox2)

		vlayout.setSpacing(3)

		groupBox.setLayout(vlayout)

		return groupBox

	def updateSigma(self):
		if self.filter_btn.isChecked():
			try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
					# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
				if self.sigma_input.text() == '':
					self.sigma = 0
				else:
					self.sigma = eval(self.sigma_input.text()) #float(self.a_input.text())
					self.sigma_input.setPlaceholderText(str(self.sigma))
					self.filter_bool = True
					self.sigma_real = self.sigma*self.L/(self.pix-1)
					self.sigma_real_label.setText("     %.2f nm"  %(self.sigma_real))
					self.plotAtoms()
			except:
				# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
				# https://www.w3schools.com/python/python_try_except.asp
				# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
				self.sigma_error = QMessageBox()
				self.sigma_error.setWindowTitle("Error")
				self.sigma_error.setText("Type in a number or numerical expression.")
				self.sigma_error.setInformativeText("Your input for sigma is: " +  str(self.sigma_input.text()))
				self.sigma_error.setIcon(QMessageBox.Warning)
				self.sigma_error.setStandardButtons(QMessageBox.Retry)
				x = self.sigma_error.exec()
		else:
			self.filter_bool = False
			self.plotAtoms()
			pass

	def initSpotifyButton(self):
		groupBox = QGroupBox(":)")
		# groupBox.setToolTip("Saves the real space image array as a .txt file,\na new .txt file with the input parameters\nand a .png of the displayed matplotlib figure")
		
		# # Add Qline Edit to enter file name
		# self.harryLabel = QLabel("HH")
		# self.harryLabel.setPixmap(QPixmap('HH.png'))
		# self.harryLabel.returnPressed.connect(self.updateSaveButton) # Connect this intput dialog whenever the enter/return button is pressed
		# self.harryLabel.setPlaceholderText('Save as...')

		# self.saveLabel = QLabel("Input file name only (no extension)")
		# Add QPushButton to open file directory whenever you want to save
		self.spotify_btn = QPushButton("Click", self)
		self.spotify_btn.clicked.connect(self.updateSpotifyButton)
		# self.spotify_btn.clicked.connect(self.openAscii)
		self.spotify_btn.setAutoDefault(False) # This is so its not default on, you can press enter on other widgets and this window wont pop up
		self.spotify_btn.setIcon(QIcon("HH.png"))
		vlayout = QVBoxLayout()
		# vlayout.addWidget(self.harryLabel)
		vlayout.addWidget(self.spotify_btn)
		groupBox.setLayout(vlayout)

		return groupBox

	def updateSpotifyButton(self):
		self.url = 'https://open.spotify.com/album/5r36AJ6VOJtp00oxSkBZ5h?si=RVB4b2eNR6iBZH4vQln6rQ'
		webbrowser.open(self.url)

	
	def ignore_this(self):
		pass
		# def initInstructionsBtn(self):
		# 	groupBox = QGroupBox()
		# 	self.instruct_btn = QPushButton("Click for instructions", self)
		# 	self.instruct_btn.clicked.connect(self.showInstructionsBtn)
		# 	self.instruct_btn.setAutoDefault(False) # This is so its not default on, you can press enter on other widgets and this window wont pop up
		# 	# self.instruct_btn.adjustSize()
		# 	# self.instruct_btn.setGeometry(100,50,500,650)
		# 	hlayout = QHBoxLayout()
		# 	hlayout.addWidget(self.instruct_btn)
		# 	groupBox.setLayout(hlayout)

		# 	return groupBox

		# def showInstructionsBtn(self, s): # Popup button syntax https://www.techwithtim.net/tutorials/pyqt5-tutorial/messageboxes/
		# 	self.instruct = QMessageBox()
		# 	self.instruct.setText("Instructions")
		# 	self.instruct.setInformativeText("Press enter after typing in an input.")
		# 	self.instruct.setDetailedText("Moire lattice: Choose yes or no. If 'No' is selected, changing the lattice 2 parameters won't affect anything until you select 'Yes'."\
		# 									"\nPixels: Creates a [pix x pix] atomic topography."	)
		# 	self.instruct.setIcon(QMessageBox.Information)
		# 	self.instruct.setStandardButtons(QMessageBox.Ok)
		# 	x = self.instruct.exec()
		# def init_vt(self):
		# 	groupBox = QGroupBox("Estimated time calculator:")
		# 	groupBox.setToolTip("Set velocity of scanner tip to estimate how long the measurement should take on the STM")
		# 	vlayout = QVBoxLayout(self)

		# 	self.calc_btn = QPushButton("Click to estimate time to take image")
		# 	self.calc_btn.clicked.connect(self.popup_calc)
		# 	self.calc_btn.setAutoDefault(False) 

		# 	vlayout.addWidget(self.calc_btn)

		# 	groupBox.setLayout(vlayout)

			

		# 	return groupBox # Adds this groupbox to the grid of widgets


		# def popup_calc(self):
		# 	self.calculator_window = CalculatorWindow(SimulatorWidget)
		# 	# self.calculator_window.setGeometry(500,150,600,600)
		# 	self.calculator_window.show()
		# 	# try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
		# 	# 		# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
		# 	# 	self.vt = eval(self.vt_input.text()) #float(self.a_input.text())
		# 	# 	self.vt_units.setText(' nm/s') 
		# 	# 	self.vt_input.setPlaceholderText(str(self.vt))
		# 	# 	# self.plotAtoms() 
		# 	# 	self.time_est = (2 * self.pix * self.L) / self.vt


		# 	# 	# self.time_mins = self.time_est / 60
		# 	# 	# self.time_hrs = self.time_est / 60 /

		# 	# 	self.time_label.setText("Estimated time: %.2f s" % (self.time_est))

		# 	# except :
		# 	# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
		# 	# https://www.w3schools.com/python/python_try_except.asp
		# 	# # Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
		# 	# self.calc_window = QMessageBox()
		# 	# # self.calc_window.setWindowTitle("Click")
		# 	# # if self.pix >= 8192:
		# 	# # self.pix_error.setText('Number of pixels exceeds 8192. Please enter a smaller number.')
		# 	# # else:
		# 	# self.calc_window.setText("Type in the tip scanner velocity.")
			
		# 	# self.vt_error.setInformativeText("Your input for vt is: " +  str(self.vt_input.text()))
		# 	# self.vt_error.setIcon(QMessageBox.Warning)
		# 	# self.vt_error.setStandardButtons(QMessageBox.Retry)

			# self.vt_label = QLabel("vt")
			# self.vt_input = QLineEdit(self)
			# self.vt_input.returnPressed.connect(self.update_vt) # Connect this intput dialog whenever the enter/return button is pressed
			# self.vt_input.setPlaceholderText(str(self.vt))

			# # Define label to display the value of the slider next to the slider
			# self.vt_units = QLabel(" nm/s", self) # Display the corrected value. only up to 2 decimal pts
			# self.vt_units.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
			# self.vt_units.setMinimumWidth(30)


			# hbox = QHBoxLayout(self)
			# hbox.addWidget(self.vt_label)
			# hbox.addWidget(self.vt_input)
			# hbox.addWidget(self.vt_units)

			# self.time_tot = (2 * self.pix * self.L) / self.vt # Total time in seconds
			# # self.time_label = QLabel("Estimated time: %.2f s" % (self.time_est), self)
			# self.hrs = self.time_tot // 3600
			# self.remain = self.time_tot - self.hrs*3600
			# self.mins = self.remain // 60
			# self.sec = self.tot - (3600*hrs) - (60*mins)

			# # total time estimated in hrs.mins.sec format
			# # print(str(hrs) + 'h ' + str(mins) + 'min ' + str(sec) + 's')

			# vlayout = QVBoxLayout(self)
			# # vlayout.addWidget(self._param_label)
			# vlayout.addLayout(hbox)
			# vlayout.addWidget(self.time_label)

			# groupBox = QGroupBox("Calculator")
			# groupBox.setLayout(vlayout)



			# x = self.calc_window.exec()


			# def update_vt(self):

			# # 	# self.vt_label = QLabel("vt")
			# # 	# self.vt_input = QLineEdit(self)
			# # 	# self.vt_input.returnPressed.connect(self.update_vt) # Connect this intput dialog whenever the enter/return button is pressed
			# # 	# self.vt_input.setPlaceholderText(str(self.vt))

			# # 	# # Define label to display the value of the slider next to the slider
			# # 	# self.vt_units = QLabel(" nm/s", self) # Display the corrected value. only up to 2 decimal pts
			# # 	# self.vt_units.setAlignment(Qt.AlignCenter | Qt.AlignVCenter)
			# # 	# self.vt_units.setMinimumWidth(30)

			# # 	# hbox = QHBoxLayout(self)
			# # 	# hbox.addWidget(self.vt_label)
			# # 	# hbox.addWidget(self.vt_input)
			# # 	# hbox.addWidget(self.vt_units)

			# # 	# self.time_est = (2 * self.pix * self.L) / self.vt
			# # 	# self.time_label = QLabel("Estimated time: %.2f s" % (self.time_est), self)


			# # 	# # vlayout.addWidget(self._param_label)
			# # 	# vlayout.addLayout(hbox)
			# # 	# vlayout.addWidget(self.time_label)

			# # 	# groupBox.setLayout(vlayout)
			# 	try: 	# I used eval() instead of float() in case an input is a mathematical expression like '3.2-1.9' 
			# 	# 		# https://stackoverflow.com/questions/9383740/what-does-pythons-eval-do
			# 		self.vt = eval(self.vt_input.text()) #float(self.a_input.text())
			# 		self.vt_units.setText(' nm/s') 
			# 		self.vt_input.setPlaceholderText(str(self.vt))
			# 	# 	# self.plotAtoms() 
			# 		self.time_est = (2 * self.pix * self.L) / self.vt


			# 	# 	# self.time_mins = self.time_est / 60
			# 	# 	# self.time_hrs = self.time_est / 60 /

			# 		self.time_label.setText("Estimated time: %.2f s" % (self.time_est))

			# 	except :
			# 	# try/except to handle errors in case the input is a string, so it doesnt just crash, instead it pops up an error window
			# 	# https://www.w3schools.com/python/python_try_except.asp
			# 	# Pop up button syntax: https://pythonprogramminglanguage.com/pyqt5-message-box/
			# 		self.vt_error = QMessageBox()
			# 		self.vt_error.setWindowTitle("Error")
			# 		# if self.pix >= 8192:
			# 		# self.pix_error.setText('Number of pixels exceeds 8192. Please enter a smaller number.')
			# 		# else:
			# 		self.vt_error.setText("Please enter a number or numerical expression.")
					
			# 		self.vt_error.setInformativeText("Your input for vt is: " +  str(self.vt_input.text()))
			# 		self.vt_error.setIcon(QMessageBox.Warning)
			# 		self.vt_error.setStandardButtons(QMessageBox.Retry)

		# def openAscii(self): # found how to open a txt file in a seperate window from here https://stackoverflow.com/questions/22829078/opening-a-text-file-in-another-window-using-push-button-in-pyqt
		# 	# openFile = QAction("Open file", self) 
		# 	# openFile.setShortcut("Ctrl+O")
		# 	# openFile.setStatusTip("Open file")
		# 	# openFile.triggered.connect(self.file_open)

		# 	open('HH_ascii.txt')


	# TO SUPPRESS QLAYOUT WARNING IN TERMINAL. from: https://stackoverflow.com/questions/25660597/hide-critical-pyqt-warning-when-clicking-a-checkboc
	def handler(msg_type, msg_log_context, msg_string):
  	  pass

	qInstallMessageHandler(handler)







