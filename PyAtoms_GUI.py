 # -*- coding: utf-8 -*-
"""
PYATOMS ATOM SIMULATOR
Created on Mon Nov 15 14:45:06 2021
@author: Asari
"""

import sys
import time

import numpy as numpy

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

# PyQT imports for creating widgets/etc
from PyQt5.QtWidgets import * #QApplication, QSplashScreen, QProgressBar, QWidget, QLabel, QPushButton, QSpinBox, QMenu, QComboBox, QMainWindow, QHBoxLayout, QVBoxLayout, QSlider, QGroupBox, QGridLayout, QRadioButton, QDialog, QLineEdit, QInputDialog
from PyQt5.QtCore import * #Qt
from PyQt5.QtGui import * #QPainter, QColor, QPixmap, QFont


# THESE TWO ARE FOR EMBEDDING MATPLOTLIB PLOTS INTO PYQT5 GUIs
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# FigureCanvasQTAgg : It will provide the canvas for the figure
# NavigationToolbar2QT : It will provide the tool bar for the matplotlib figs (zooming in, panning, saving image, etc)
# https://www.geeksforgeeks.org/how-to-embed-matplotlib-graph-in-pyqt5/




# IMPORT OTHER FILES:
from PyAtoms_Widgets import SimulatorWidget
from hexatoms import hexatoms
from squareatoms import squareatoms
from moirelattice import moirelattice


# main window
# which inherits QDialog
class Window(QDialog):
       
    # constructor - always include this
    def __init__(self, parent=None):
        super(Window, self).__init__(parent) # always need to do this because our class Window is inheriting from parent QDialog....
        self.initUI() # Call the function that initializes everything
        
        

    def initUI(self):
        # Define default instance variables 
        # Use 'self.' when defining basically any thing in the class - these are class member variables (?) they can be accessed in all functions throughout the class

        # Call functions to initialize everything else
        self.setWindowTitle('PyAtoms v. 1.0') # Sets the title on the external window that pops up when you run the code
        self.initGeo() # Sets size of the popup gui window
        self.initWidgetsGrid() # For placing multiple widgets in the popup gui in a grid layout
        self.show()



    ## Tried to change the numbers here to change the size of the window but it doesnt work?
    def initGeo(self): 
        # Set geometry of popup gui window
        self.setGeometry(100, 100, 1250,500)#self.width(),self.height())#1200, 850)
        # self.setStyleSheet("background: gray;") # Change color of background in window
        # self.setStyleSheet("color:red") # Change color of allll the displayed text in the gui
        # self.setStyleSheet("color: magenta;background: gray")
        # self.setStyleSheet("border: 1px solid black;")

      
        # Hard-coded these numbers just by running the code and seeing how it looked when I changed the numbers
        # These numbers set the size of the window that pops up
        self.x = self.width() // 3 + 20 
        self.y = 40
        self.w = (2*self.width()) // 3 - 100
        self.h = self.height() - 100
        self.SimWidget = SimulatorWidget(self, self.x,self.y,self.w,self.h) # Define an instance of the SimulatorWidget class, which is being imported thus can be accessed in this file
        self.showMaximized() # To open the window fully maximized https://www.geeksforgeeks.org/pyqt5-how-to-open-window-in-maximized-format/


    # Create a layout to place all the widgets/groupboxes in a grid layout
    # This is basically the function that calls all the other init_something functions in the child class to place them as widgets in the gui grid
    def initWidgetsGrid(self):
        # SYNTAX for adding widget to gui:
        #(addWidget(QWidget, int r, int c, int rowspan, int columnspan)) ---  Adds widget at specified row and column and having specified width and/or height
        #grid.addWidget(row, column, width, height)
        
        grid = QGridLayout()
        grid.setSpacing(11)

        # Settings for all lattices
        grid.addWidget(self.SimWidget.initMoireBtn(), 0, 1)
        grid.addWidget(self.SimWidget.initColormapDropdown(),1,1)
        # grid.addWidget(self.SimWidget.initSpotifyButton(),9, 0, 1,1)

        grid.addWidget(self.SimWidget.initImageParameters(), 0,2,2,1)

        
        grid.addWidget(self.SimWidget.initFiltering(), 0, 3)
        grid.addWidget(self.SimWidget.initSaveButton(),1, 3)

        # Add calculator widgets
        grid.addWidget(self.SimWidget.initCalcWidget(),0,4)
        grid.addWidget(self.SimWidget.initMapCalcWidget(), 1, 4)

        # Settings for first lattice
        grid.addWidget(self.SimWidget.initLattice1Parameters(), 0, 0, 3, 1)


        # Settings for second lattice
        grid.addWidget(self.SimWidget.initLattice2Parameters(), 3, 0,3,1)

        # Lattice 3 parameters
        grid.addWidget(self.SimWidget.initLattice3Parameters(), 6, 0,3,1)


        # Add matplotlib fig/toolbar to gui
        grid.addWidget(self.SimWidget.initMatplotlibFig(),2,1,8,5) # Make the matplotlib canvas/figure the largest widget
        
        
        self.setLayout(grid)
        

    # Overriding keyPressEvent so that if the escape button is pressed, it doesn't automatically close the program
    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            # print("Esc clicked")
            self.esc_error = QMessageBox()
            self.esc_error.setWindowTitle("Warning!")
            self.esc_error.setText("Escape button pressed. Do you want to exit the program?")
            self.esc_error.setIcon(QMessageBox.Warning)
            self.esc_error.setStandardButtons(QMessageBox.Cancel)
            self.esc_error.addButton(QPushButton("Exit program"), QMessageBox.NoRole) # Add a custom button, from https://stackoverflow.com/questions/15682665/how-to-add-custom-button-to-a-qmessagebox-in-pyqt4
            x = self.esc_error.exec()
            # print((x))
            if x == 0: # if x =0, this means the Exit program button was clicked (found out by just printing the value of x for each button)
                sys.exit(app.exec_()) # Exit program if the exit button was clicked




# Handle high resolution displays:  https://stackoverflow.com/questions/39247342/pyqt-gui-size-on-high-resolution-screens
if hasattr(Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    
# To run, go to the file path location in terminal and type  'python PyAtoms_GUI.py'
# driver code
if __name__ == '__main__':  # this won't be run when imported... https://stackoverflow.com/questions/6523791/why-is-python-running-my-module-when-i-import-it-and-how-do-i-stop-it
                            # # code here will only run when you invoke 'python main.py'
    # creating apyqt5 application
    app = QApplication(sys.argv)

    # Create splash screen (loading screen) from https://gist.github.com/345161974/8897f9230006d51803c987122b3d4f17
    # splash_pix = QPixmap("HH.png")
    splash_pix = QPixmap("logo_magma_Small.png")
    # splash_pix.scaledToHeight(240, Qt.SmoothTransformation)
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint) # WindowStaysOnTopHint: to keep it above all the other windows on the desktop.
    # splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint) # I commented this out and it fixed the issues of it reappearing after clicking away and back from hte gui
    splash.setEnabled(False)

    fontsize=40
    font = splash.font()
    font.setPixelSize(fontsize)
    font.setWeight(QFont.Bold)
    splash.setFont(font)

    # adding progress bar
    progressBar = QProgressBar(splash)
    progressBar.setMaximum(10)
    progressBar.setGeometry(0, splash_pix.height() - 50, splash_pix.width(), 20)

    splash.show()
    splash.showMessage("<h1><font color='white'>Loading...</font></h1>",Qt.AlignBottom, Qt.black)
    

    for i in range(1, 12):
        progressBar.setValue(i)
        t = time.time()
        while time.time() < t + 0.1:
           app.processEvents()

    # Simulate something that takes time
    time.sleep(1)

    # creating a window object
    main = Window()
       
    # showing the window
    #    main.show()

    splash.finish(main) # Remove the splash when the_editor has finished setting itself up.

    splash.close()

    # loop
    sys.exit(app.exec_())

    
    




































