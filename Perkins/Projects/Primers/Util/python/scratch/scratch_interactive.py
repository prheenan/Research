"""
Demonstrates some customized mouse interaction by drawing a crosshair that follows 
the mouse.


"""

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
from pyqtgraph.Point import Point

#generate layout
app = QtGui.QApplication([])
win = pg.GraphicsWindow()
win.setWindowTitle('pyqtgraph example: crosshair')
label = pg.LabelItem(justify='right')
win.addItem(label)
p1 = win.addPlot(row=1, col=0)
p2 = win.addPlot(row=2, col=0)

region = pg.LinearRegionItem()
region.setZValue(10)
# Add the LinearRegionItem to the ViewBox, but tell the ViewBox to exclude this 
# item when doing auto-range calculations.
p2.addItem(region, ignoreBounds=True)

#pg.dbg()
p1.setAutoVisible(y=True)


#create numpy arrays
#make the numbers large to show that the xrange shows data from 10000 to all the way 0
n = int(1e6)
data1 = 10000 + 15000 * pg.gaussianFilter(np.random.random(size=n), 10) + 3000 * np.random.random(size=n)
data2 = 15000 + 15000 * pg.gaussianFilter(np.random.random(size=n), 10) + 3000 * np.random.random(size=n)

p1.plot(data1, pen="r")
p1.plot(data2, pen="g")

p2.plot(data1, pen="w")

def update():
    region.setZValue(10)
    minX, maxX = region.getRegion()
    p1.setXRange(minX, maxX, padding=0)    

region.sigRegionChanged.connect(update)

def updateRegion(window, viewRange):
    rgn = viewRange[0]
    region.setRegion(rgn)

p1.sigRangeChanged.connect(updateRegion)

region.setRegion([1000, 2000])

#cross hair
vLine = pg.InfiniteLine(angle=90, movable=False)
hLine = pg.InfiniteLine(angle=0, movable=False)
p1.addItem(vLine, ignoreBounds=True)
p1.addItem(hLine, ignoreBounds=True)


vb = p1.vb

# see: https://groups.google.com/forum/#!msg/pyqtgraph/nI5hZJ85N7M/TsPoJfk0V8kJ
# prints all PyQtGraph Signals
def PrintPyQtGraphSignals():
    import pyqtgraph as pg
    import inspect
    for n,o in pg.__dict__.items():
        if inspect.isclass(o) and issubclass(o, pg.QtCore.QObject):
            for m,p in o.__dict__.items():
                if 'unbound signal' in str(p):
                    print n, m 

def mouseMoved(evt):
    pos = evt[0]  ## using signal proxy turns original arguments into a tuple
    if p1.sceneBoundingRect().contains(pos):
        mousePoint = vb.mapSceneToView(pos)
        x = mousePoint.x()
        y = mousePoint.y()
        index = int(x)
        if index > 0 and index < len(data1):
            label.setText("<span style='font-size: 12pt'>x=%0.1f,   <span style='color: red'>y1=%0.1f</span>,   <span style='color: green'>y2=%0.1f</span>" % (mousePoint.x(), data1[index], data2[index]))
        vLine.setPos(x)
        hLine.setPos(y)
        print(x)

def mouseClick(evt):
    print("Click!")
    print(evt)

PrintPyQtGraphSignals()
updateHz=30
proxy = pg.SignalProxy(p1.scene().sigMouseClicked, rateLimit=updateHz,
                       slot=mouseClick)
proxy = pg.SignalProxy(p1.scene().sigMouseMoved, rateLimit=updateHz,
                       slot=mouseMoved)
#p1.scene().sigMouseMoved.connect(mouseMoved)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
