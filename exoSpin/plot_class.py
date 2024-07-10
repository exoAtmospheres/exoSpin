'''
ExoSpin run script - Plot Class


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''



# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports
import sys, os, glob
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as c
import pickle

# ------------------------------------------------------------------------------------------------------------------------------------------------

class Plot ():

    def __init__(self , type_ , x , y , xlabel , ylabel , color , title):

        self.type = type_
        self.x = x
        self.y = y
        self.xlabel= xlabel
        self.ylabel = ylabel
        self.color = color
        self.title = title
        self.bins = 200

    def set_color(self,new_color):
        self.color = new_color

    def set_title(self,new_title):
        self.title = new_title

    def set_type(self,new_color):
        self.color = new_color

    def set_x(self,new_x):
        if new_x.size == self.x.size:
            self.x = new_x
        else:
            return 'Error of size : Your x-vector hasn\'t the good size'

    def set_x(self,new_y):
        if new_y.size == self.y.size:
            self.y = new_y
        else:
            return 'Error of size : Your y-vector hasn\'t the good size'
    
    def set_xlabel(self, label):
        self.xlabel = label
    
    def set_ylabel(self, label):
        self.ylabel = label

    def plot(self):
        fig = plt.figure()
        if self.type == 'Histogram':
            y, x, _ = plt.hist(self.x, bins=self.bins, density=True, color=self.color)
            plt.title(self.title)
            plt.xlabel(self.xlabel)
            plt.show()
        else : 
            fig = plt.plot(self.x,self.y,color=self.color)
            plt.title(self.title)
            plt.xlabel(self.xlabel)
            plt.ylabel(self.ylabel)
            plt.legend()
            plt.show()
        