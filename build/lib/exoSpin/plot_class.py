'''
ExoSpin - Plot Class


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
        """
        Initialize every parameters for the Plot class object.

        Args:
            type_ (String) : a String representing the type of plot.
                            {
                                'Histogram'
                                'PDF'
                            }
            x (numpy.ndarray) : 1D array representing x-values of the plot.
            y (numpy.ndarray) : 1D array representing y-values of the plot.
            xlabel (String) : a String representing the label of the x-axis.
            ylabel (String) : a String representing the label of the y-axis.
            color (String) : a String representing the color of the plot.
            title (String) : a String representing the titlte of the plot.

        Attributes:
            type_ (String) : a String representing the type of plot.
                            {
                                'Histogram'
                                'PDF'
                            }
            x (numpy.ndarray) : 1D array representing x-values of the plot.
            y (numpy.ndarray) : 1D array representing y-values of the plot.
            xlabel (String) : a String representing the label of the x-axis.
            ylabel (String) : a String representing the label of the y-axis.
            color (String) : a String representing the color of the plot.
            title (String) : a String representing the titlte of the plot.
            bins (int) : a integer that sets the histogram bins.

        """

        self.type = type_
        self.x = x
        self.y = y
        self.xlabel= xlabel
        self.ylabel = ylabel
        self.color = color
        self.title = title
        self.bins = 200

    def set_color(self,new_color):
        """
        Set a new color for the plot.

        Args:
            new_color (String): a String that represents the new color for the plot.

        Raises:
            ValueError: If input and self attribute don't have the same type.
        """
        if type(self.color)!=type(new_color):
            raise ValueError("The input and self attribute don't have the same type")

        self.color = new_color

    def set_title(self,new_title):
        """
        Set a new title for the plot.

        Args:
            new_color (String): a String that represents the new title for the plot.
            
        Raises:
            ValueError: If input and self attribute don't have the same type.
        """

        if type(self.title)!=type(new_title):
            raise ValueError("The input and self attribute don't have the same type")

        self.title = new_title

    def set_x(self,new_x):
        """
        Set a new x-axis for the plot

        Args:
            new_x (numpy.ndarray): 1D array that represents the new x-axis for the plot.
            
        Raises:
            ValueError: If input and self attribute don't have the same type and length.
        """

        if type(self.x)!=type(new_x):
            raise ValueError("The input and self attribute don't have the same type")
        if self.x.size != new_x.size:
            raise ValueError("The input and self attribute don't have the same type")
        self.x = new_x

    def set_y(self,new_y):
        """
        Set a new y-axis for the plot

        Args:
            new_y (numpy.ndarray): 1D array that represents the new y-axis for the plot.
            
        Raises:
            ValueError: If input and self attribute don't have the same type and length.
        """

        if type(self.y)!=type(new_y):
            raise ValueError("The input and self attribute don't have the same type")
        if self.y.size != new_y.size:
            raise ValueError("The input and self attribute don't have the same type")
        self.y = new_y
    
    def set_xlabel(self, label):
        """
        Set a new x-axis label for the plot.

        Args:
            label (String): a String that represents the new label for the plot.
            
        Raises:
            ValueError: If input and self attribute don't have the same type.
        """

        if type(self.xlabel)!=type(label):
            raise ValueError("The input and self attribute don't have the same type")

        self.xlabel = label
    
    def set_ylabel(self, label):
        """
        Set a new y-axis label for the plot.

        Args:
            label (String): a String that represents the new label for the plot.
            
        Raises:
            ValueError: If input and self attribute don't have the same type.
        """

        if type(self.ylabel)!=type(label):
            raise ValueError("The input and self attribute don't have the same type")

        self.ylabel = label

    def plot(self):
        """
        Show the plot.
        """

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
            plt.show()
        