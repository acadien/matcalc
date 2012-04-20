#!/usr/bin/python

#Some of this code taken from:
#http://www.uppmax.uu.se/Members/andersh/vtk-with-python/isosurfaces/advanced-isosurfaces
#All of that code is GNU-GPLed... so you need to dump this code and start from scratch.

import sys,operator,os
from numpy import *
import scipy
from vtk import *
#mine
from chgcarIO import readchgcar

def usage():
    print "Usage:"
    print "%s <chg-vtk-file>"%sys.argv[0]
    sys.exit()

if not(len(sys.argv) in [2]):
    usage()

vtkfile=sys.argv[1]
#---------------------------------
#Below is taken from website
#---------------------------------
# Interaction with an Isosurface visualization

# A reader
reader = vtkStructuredPointsReader()
reader.SetFileName(vtkfile)

renWin = vtkRenderWindow()
renWin.SetWindowName("%s"%vtkfile)
renWin.SetPolygonSmoothing(1)
#renWin.SetAAFrames(2)

renWin.SetSize(500, 500)

# Create an outline of the dataset (Cube lines)
outline = vtkOutlineFilter()
outline.SetInputConnection(reader.GetOutputPort())

outlineMapper = vtkPolyDataMapper()
outlineMapper.SetInputConnection(outline.GetOutputPort())

outlineActor = vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0.0,0.0,0.5)

# Color lookup table
lut=vtkColorTransferFunction()
lut.AddRGBPoint(0.0,0,0,0)
lut.AddRGBPoint(0.850,1,0,1)
lut.AddRGBPoint(0.890,0,0,1)
lut.AddRGBPoint(0.925,0,1,0)
lut.AddRGBPoint(1.000,1,1,0)
lut.AddRGBPoint(5.0,1,1,1)

# Define initial iso value 
isovalue = 0.9

# The contour filter
isosurface = vtkContourFilter()

isosurface.SetInputConnection(reader.GetOutputPort())
isosurface.SetValue(0,isovalue)

isosurfaceMapper = vtkPolyDataMapper()
isosurfaceMapper.SetLookupTable(lut)
isosurfaceMapper.SetInputConnection(isosurface.GetOutputPort())

isosurfaceActor = vtkActor()
isosurfaceActor.SetMapper(isosurfaceMapper)

# Renderer and render window 
ren = vtkRenderer()
ren.SetBackground(1.,1.,1.)

#Create some text (the iso value)
textActor = vtkTextActor()
tp = vtkTextProperty()
tp.BoldOn()
tp.ShadowOn()
tp.ItalicOn()
tp.SetColor(1.0,0.2,0.3)
tp.SetFontFamilyToArial()
tp.SetFontSize(20)
textActor.SetTextProperty(tp)
textActor.SetInput(str(isovalue))

# Add the actors
ren.AddActor(outlineActor)
ren.AddActor(isosurfaceActor)
ren.AddActor(textActor)


renWin.AddRenderer(ren)

# Python function for the keyboard interface
def Keypress(obj, event):
    global isovalue, renWin
    key = obj.GetKeySym()
    if key == "x":
        print outlineActor.GetMapper().GetInputConnection()
    #Alter the iso-surface
    if key == "m":
        isovalue = isovalue - 0.01
    elif key == "comma":
        isovalue = isovalue - 0.001
    elif key == "period":
        isovalue = isovalue + 0.001
    elif key == "slash":
        isovalue = isovalue + 0.01
    if key in ["m","comma","period","slash"]:
        isosurface.SetValue(0,isovalue)
        textActor.SetInput("%4.3f" %(isovalue))
        tp.SetColor(lut.GetColor(isovalue))
        renWin.Render()
    
#List the keyboard shortcuts
print "Keyboard Shortcuts:"
print "m : isosurface - big step"
print ", : isosurface - small step"
print ". : isosurface + small step"
print "/ : isosurface + big step"

# add keyboard interface, initialize, and start the interactor

# Render window interactor
iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

tpc = textActor.GetPositionCoordinate()
tpc.SetCoordinateSystemToNormalizedViewport()
tpc.SetValue(0.75,0.9)

iren.AddObserver("KeyPressEvent", Keypress)
iren.Initialize()
iren.Start()



