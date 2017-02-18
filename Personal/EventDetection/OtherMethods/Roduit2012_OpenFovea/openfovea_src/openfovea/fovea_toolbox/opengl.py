#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
    This opengl code was higly inspired from the pyopengl documentation, and
    particularly on the NeHe tutorials, thank's a lot to them !
"""
from copy import deepcopy

import numpy
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from gtk.gtkgl.apputils import *

from matplotlib.pylab import pcolor
import time

from misc import ismasked

class Topography(GLScene,
             GLSceneButton,
             GLSceneButtonMotion):
    """
        Create a 3D view of an array, the color is based on the color array and
        the topography on the array.
        
        To create a view, follow the following simple steps :
        
        >>> import numpy
        >>> array = numpy.random.randn(32, 32)
        >>> color_array = numpy.random.randn(32, 32)
        >>> Topography(array, color_array)
    """
    
    def __init__(self, array, color_array=None, clut=None):
        GLScene.__init__(self,
                         gtk.gdkgl.MODE_RGB   |
                         gtk.gdkgl.MODE_DEPTH |
                         gtk.gdkgl.MODE_DOUBLE)
        self.window = 0 # the glut window
        self.rotation = 0.0 # the rotation angle for the object
        self.array = array
        self.rot = [1, 1, 1]
        self.zoom = 20
        self.translate = [0, 0, 0]
        # This is the Color Look-Up Table (CLUT)
        if color_array is not None:
            if clut is None:
                clut = (array.min(), array.max())
            self.clut = pcolor(numpy.array([[clut[0], clut[1]] for i in range(2)]))
        else:
            self.clut = None # Have to define it later...
        
        self.mouse_eventpos_begin = [0, 0] # To record the position 
                                           # of the mouse during an event
        self.slice = {'begin' : [0, 0, 0],  # This is the portion of the array
                      'end' : [-1, -1, -1]} # to display
        ### Initialize the Window
        # pass argument to init
        if array is not None:
            self.coordinate = make_coordinate(array, self.clut)
            self.all_coordinate = deepcopy(self.coordinate)
        glHint(GL_CLIP_VOLUME_CLIPPING_HINT_EXT, GL_FASTEST)
        
        glutInit()
        #glutCreateWindow('OpenGL Tutorial 1')
    def init(self):
        glEnable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)
        
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.5, 0.5, 0.8))
        glLightfv(GL_LIGHT0, GL_POSITION, (0, 2, 2, 0))
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
    def keyboard_perss(self, event):
        print event
    def button_press(self, width, height, event):
        if event.state & gtk.gdk.BUTTON4_MASK:
            print "ZOOM"
        if event.state & gtk.gdk.BUTTON5_MASK:
            print "UNZOOM"
        if event.state & gtk.gdk.SCROLL_MASK:
            print "Scroll"
        self.mouse_eventpos_begin = [event.x, event.y]
    
    def button_release(self, width, height, event):
        pass
    def button_motion(self, width, height, event):
        if event.state & gtk.gdk.BUTTON1_MASK:
            deg_factor = numpy.pi/180
            rot_x = ((event.x - self.mouse_eventpos_begin[0]) / width) * 360
            rot_y = ((event.y - self.mouse_eventpos_begin[1]) / height) * 360
            self.rot[0] = (self.rot[0] + rot_y) % 360
            self.rot[1] = (self.rot[1] + numpy.cos(self.rot[0] * deg_factor) * rot_x) % 360
            self.rot[2] = (self.rot[2] - numpy.sin(self.rot[0] * deg_factor) * rot_x) % 360
        if event.state & gtk.gdk.BUTTON2_MASK:
            self.translate[0] = self.translate[0] + (
                                event.x - self.mouse_eventpos_begin[0]) / width
            self.translate[1] = self.translate[1] - (
                                event.y - self.mouse_eventpos_begin[1]) / height
        self.mouse_eventpos_begin = [event.x, event.y]
        self.invalidate()
    def rotate(self):
        """
            The rotation of the array
        """
        
        self.rotx += 2
        if self.rotx > 360:
            self.rotx = self.rotx - 360
        self.roty += 2
        if self.roty > 360:
            self.roty = self.roty - 360
        glutPostRedisplay()
            
    def reshape(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if width > height:
            w = float(width) / float(height)
            glFrustum(-w, w, -1.0, 1.0, 5.0, 60.0)
        else:
            h = float(height) / float(width)
            glFrustum(-1.0, 1.0, -h, h, 5.0, 60.0)
        #glMatrixMode(GL_MODELVIEW)
    def display(self, width, height):
        """
            The main drawing funtion
        """
        # Clear The Screen And The Depth Buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()    # Reset The View
        glMatrixMode(GL_PROJECTION)
        # The text
        glColor3f(1.0, 1, 1)
        glRasterPos2f(-0.9,0.9)
        to_display = "x : [%i : %i], y : [%i : %i], z : [%i : %i]"%(
                     self.slice['begin'][0], self.slice['end'][0],
                     self.slice['begin'][1], self.slice['end'][1],
                     self.slice['begin'][2], self.slice['end'][2],
                     )
        for string in to_display:
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13, ord(string))
        # Modify the projection matrix
        gluPerspective(self.zoom, 1., 1., 500.0)
        glTranslatef(0.0,  0.0, -4)
        glRotate(self.rot[0],1.0,0.0,0.0)
        glRotate(self.rot[1], 0.0, 1.0, 0.0)
        glRotate(self.rot[2], 0.0, 0.0, 1.0)
        glTranslate(self.translate[0], self.translate[1], self.translate[2])
        # Every thing after is modified by the projection matrix
        self.generate_vertex(self.coordinate)
            
        glFlush()
        
    def generate_vertex(self, coordinate):
        glBegin(GL_QUADS)
        for item in coordinate:
            glNormal3f(item['normal'][0], item['normal'][1], item['normal'][2])
            # Bottom left
            glColor3f(item['color'][0][0],
                      item['color'][0][1],
                      item['color'][0][2])
            glVertex3f(item['vertex_coord'][0][0],
                       item['vertex_coord'][0][1],
                       item['vertex_coord'][0][2])
            # upper left
            glColor3f(item['color'][1][0],
                      item['color'][1][1],
                      item['color'][1][2])
            glVertex3f(item['vertex_coord'][1][0], 
                       item['vertex_coord'][1][1],
                       item['vertex_coord'][1][2])
            # bottom right
            glColor3f(item['color'][3][0],
                      item['color'][3][1],
                      item['color'][3][2])
            glVertex3f(item['vertex_coord'][3][0],
                       item['vertex_coord'][3][1],
                       item['vertex_coord'][3][2])
            # upper right
            glColor3f(item['color'][2][0],
                      item['color'][2][1],
                      item['color'][2][2])
            glVertex3f(item['vertex_coord'][2][0],
                       item['vertex_coord'][2][1],
                       item['vertex_coord'][2][2])
        glEnd()
        glLineWidth(1.5)
        glBegin(GL_LINES)
        for item in coordinate:
            glColor3f(0., 0., 0.)
            glVertex3f(item['vertex_coord'][0][0],
                       item['vertex_coord'][0][1],
                       item['vertex_coord'][0][2])
            glVertex3f(item['vertex_coord'][1][0],
                       item['vertex_coord'][1][1],
                       item['vertex_coord'][1][2])
            glColor3f(0., 0., 0.)
            glVertex3f(item['vertex_coord'][1][0],
                       item['vertex_coord'][1][1],
                       item['vertex_coord'][1][2])
            glVertex3f(item['vertex_coord'][3][0],
                       item['vertex_coord'][3][1],
                       item['vertex_coord'][3][2])
            glColor3f(0., 0., 0.)
            glVertex3f(item['vertex_coord'][2][0],
                       item['vertex_coord'][2][1],
                       item['vertex_coord'][2][2])
            glVertex3f(item['vertex_coord'][0][0],
                       item['vertex_coord'][0][1],
                       item['vertex_coord'][0][2])
            glColor3f(0., 0., 0.)
            glVertex3f(item['vertex_coord'][3][0],
                       item['vertex_coord'][3][1],
                       item['vertex_coord'][3][2])
            glVertex3f(item['vertex_coord'][2][0],
                       item['vertex_coord'][2][1],
                       item['vertex_coord'][2][2])
        glEnd()
class Tomography(Topography):
    def __init__(self, array, clut):
        Topography.__init__(self, None)
        self.array = array
        if clut is None:
            clut = (array.min(), array.max())
        clut = numpy.array([[clut[0], clut[1]] for i in range(2)])
        self.clut = pcolor(clut)
        self.all_coordinate = make_tomo_coordinate(array, self.clut)
        self._slice_it()
    def slice_it(self, dirct, qty):
        """
            Slice the array in a certain direction in a certain quantity.
            
            parameters :
                dirct : str
                    The direction to slice.
                    'x' to slice in x from the front.
                    'X' to slice in x from the back.
                    'y' to slice in y from the left.
                    'Y' to slice in y from the right.
                    'z' to slice in z from the top.
                    'Z' to slice in z from the bottom.
        """
        if dirct == 'x' and (0 <= self.slice['begin'][0] + 
                             qty <= self.array.shape[0]):
            self.slice['begin'][0] += qty
        elif dirct == 'X' and (0 <= self.slice['end'][0] + 
                               qty <= self.array.shape[0]):
            self.slice['end'][0] += qty
        elif dirct == 'y' and (0 <= self.slice['begin'][1] + 
                               qty <= self.array.shape[1]):
            self.slice['begin'][1] += qty
        elif dirct == 'Y' and (0 <= self.slice['end'][1] + 
                               qty <= self.array.shape[1]):
            self.slice['end'][1] += qty
        elif dirct == 'z' and (0 <= self.slice['begin'][2] + 
                               qty <= self.array.shape[2]):
            self.slice['begin'][2] += qty
        elif dirct == 'Z' and (0 <= self.slice['end'][2] + 
                               qty <= self.array.shape[2]):
            self.slice['end'][2] += qty
        self._slice_it()
    def _slice_it(self):
        """
            Slice the array.
            
            For example : self.slice_it([0, 0, 4], [-1, -1, -1])
            will slice the array in z by removing the 4 (0 to 3) first elements.
            
            Parameters :
                begin : array_like
                    begin[0] is the first x
                    begin[1] the first y
                    begin[2] the first z
                end : array_like
                    end[0] is the last x
                    end[1] is the last y
                    end[2] is the last z
            
            Retruns :
                nothing. Just the self.coordinate are modified according to it.
        """
        
        factor = float(max(self.array.shape))
        _begin = [0, 0, 0]
        _end = [-1, -1, -1]
        for i in range(3):
            if self.slice['end'][i] == -1:
                self.slice['end'][i] = self.array.shape[i]
            _begin[i] = ((float(self.slice['begin'][i]) / factor) - 
                          self.array.shape[i]/(2*factor))
            _end[i] = ((float(self.slice['end'][i]) / factor) - 
                        self.array.shape[i]/(2*factor))
        self.coordinate = list()
        for item in self.all_coordinate:
            # find items that are outside the slice.
            eliminate = 0
            for vertex in item['vertex_coord']:
                if not eliminate and (
                                vertex[0] < _begin[0] or vertex[0] > _end[0] or
                                vertex[1] < _begin[1] or vertex[1] > _end[1] or
                                vertex[2] < _begin[2] or vertex[2] > _end[2]):
                    eliminate = 1
            if not eliminate:
                self.coordinate.append(item)
        # Add the borders...
        borders = generate_array_borders(self.array,
                                         self.slice['begin'], self.slice['end'],
                                         self.clut)
        self.coordinate = self.coordinate + \
                          [i for i in borders if i is not None]
def make_tomo_coordinate(array, color_map):
    """
        This function creates the coordinate of the vertex and the corresponding
        color.
        The array is a 3D array which represent the space with pixels coordinate
        corresponding to the location in the array and the color value coded in
        that location.
        
        Parameters :
            array : numpy.ma.array
                This is an array with mask. Masked values are transparent
                pixels, otherwise, the pixel values represent the color.
        
        Returns :
            coordinate : list
                List of vertex to display.
    """
    #array = array/array.max() # to be between 0 and 1.
    #plot_array = numpy.array([[array.min(), array.max()] for i in range(2)])
    #plot_array = pcolor(plot_array) # to have a cmap
    coordinate = list()
    #non_zero_indice = array.mask.nonzero()
    for px, y, z in zip(array.mask.nonzero()[0],
                       array.mask.nonzero()[1],
                       array.mask.nonzero()[2]):
        item = generate_border_vertex(array, px, y, z, color_map)
        coordinate += [i for i in item if i is not None]
                
    # generage the top face and bottom face:
    [x_list, y_list] = numpy.meshgrid(range(array.shape[0]),
                                      range(array.shape[1]))
    for px, y in zip(x_list.flat, y_list.flat):
        item = [generate_vertex(array, px, y, array.shape[2]-1, color_map),
               generate_vertex(array, px, y, 0, color_map)]
        coordinate += [i for i in item if i is not None]
    # generate the front and back faces:
    coord_list = numpy.meshgrid(range(array.shape[0]), range(array.shape[2]))
    for px, z in zip(coord_list[0].flat, coord_list[1].flat):
        item = [generate_vertex(array, px, array.shape[1]-1, z,
                                color_map, face = 'front'),
                generate_vertex(array, px, 0, z,
                                color_map, face = 'front')]
        coordinate += [i for i in item if i is not None]
    # generate the right and left faces:
    coord_list = numpy.meshgrid(range(array.shape[1]), range(array.shape[2]))
    for y, z in zip(coord_list[0].flat, coord_list[1].flat):
        item = [generate_vertex(array, array.shape[0]-1, y, z,
                                color_map, face = 'right') ,
                generate_vertex(array, 0, y, z,
                                color_map, face = 'right')]
        coordinate += [i for i in item if i is not None]
    return coordinate
def generate_array_borders(array, begin, end, color_map):
    """
        Generate the 
    """
    #array = array/array.max() # to be between 0 and 1.
    #plot_array = numpy.array([[array.min(), array.max()] for i in range(2)])
    #plot_array = pcolor(plot_array)
    coordinate = list()
    if begin[0] != 0:
        # regenerate the front face
        coord_list = numpy.meshgrid(xrange(begin[1], end[1]),
                                    xrange(begin[2], end[2]))
        coordinate += [generate_vertex(array, begin[0], y, z,
                                       color_map, face = 'right')
                        for y, z in zip(coord_list[0].flat, coord_list[1].flat)]
    if end[0] != array.shape[0]:
        # regenerate the back face
        coord_list = numpy.meshgrid(xrange(begin[1], end[1]),
                                    xrange(begin[2], end[2]))
        coordinate += [generate_vertex(array, end[0], y, z,
                                       color_map, face = 'right')
                        for y, z in zip(coord_list[0].flat, coord_list[1].flat)]
    if begin[1] != 0:
        # regenerate the left face
        coord_list = numpy.meshgrid(xrange(begin[0], end[0]),
                                    xrange(begin[2], end[2]))
        coordinate += [generate_vertex(array, px, begin[1], z,
                                       color_map, face = 'front')
                        for px, z in zip(coord_list[0].flat, coord_list[1].flat)]
    if end[1] != array.shape[1]:
        # regenerate the left face
        coord_list = numpy.meshgrid(xrange(begin[0], end[0]),
                                    xrange(begin[2], end[2]))
        coordinate += [generate_vertex(array, px, end[1], z,
                                       color_map, face = 'front')
                        for px, z in zip(coord_list[0].flat, coord_list[1].flat)]
    if begin[2] != 0:
        # regenerate the left bottom
        coord_list = numpy.meshgrid(xrange(begin[0], end[0]),
                                    xrange(begin[1], end[1]))
        coordinate += [generate_vertex(array, px, y, begin[2],
                                       color_map, face = 'floor')
                        for px, y in zip(coord_list[0].flat, coord_list[1].flat)]
    if end[2] != array.shape[2]:
        # regenerate the left top
        coord_list = numpy.meshgrid(xrange(begin[0], end[0]),
                                    xrange(begin[1], end[1]))
        coordinate += [generate_vertex(array, px, y, end[2],
                                       color_map, face = 'floor')
                        for px, y in zip(coord_list[0].flat, coord_list[1].flat)]
    return coordinate
def generate_border_vertex(array, px, y, z, color_map):
    """
        Generates the vertex that touch this one.
    """
    # The vertex is a 3d vertex composed on 6 faces
    #
    #       (x, y+1, z+1)_______
    #                  /·     /| (x+1, y+1, z+1)
    #                 / ·    / |
    #    (x, y, z+1) +------+  |
    #                |  ·······+ (x+1, y+1, z)
    #   (x, y+1, z)  | ·    | /      
    #                |·     |/       
    #                +------+ (x+1, y, z)
    #                (x, y, z)
    block = {'right' : px + 1 == array.shape[0],
             'left' : px == 0,
             'front' : y + 1 == array.shape[1],
             'back' : y == 0,
             'top' : z + 1 == array.shape[2],
             'bottom' : z == 0}
    pixels = list()
    factor = float(max(array.shape))
    pos = [((float(px) / factor) - array.shape[0] / (2 * factor)),
           ((float(y) / factor) - array.shape[1] / (2 * factor)),
           ((float(z) / factor) - array.shape[2] / (2 * factor))]
    pos_u = [((float(px + 1) / factor) - array.shape[0] / (2 * factor)),
           ((float(y+1) / factor) - array.shape[1] / (2 * factor)),
           ((float(z+1) / factor) - array.shape[2] / (2 * factor))]
    pos_d = [((float(px - 1) / factor) - array.shape[0] / (2 * factor)),
           ((float(y-1) / factor) - array.shape[1] / (2 * factor)),
           ((float(z-1) / factor) - array.shape[2] / (2 * factor))]
    if not block['bottom'] and not array.mask[px, y, z-1]:
        # Bottom pixel exist... create the bottom face.
        color = color_map.to_rgba(array[px, y, z-1])
        pixels.append({'vertex_coord' : [[pos[0], pos[1], pos[2]],
                                              [pos[0], pos_u[1], pos[2]],
                                              [pos_u[0], pos[1], pos[2]],
                                              [pos_u[0], pos_u[1], pos[2]]],
                            'color' : [color, color, color, color],
                            'normal' : NORMAL['bottom']})
                            
    if not block['right'] and not array.mask[px + 1, y, z]:
        # Right pixel exist... create the right face.
        color = color_map.to_rgba(array[px + 1, y, z])
        pixels.append({'vertex_coord' : [[pos_u[0], pos[1], pos[2]],
                                        [pos_u[0], pos_u[1], pos[2]],
                                        [pos_u[0], pos[1], pos_u[2]],
                                        [pos_u[0], pos_u[1], pos_u[2]]],
                        'color' : [color, color, color, color],
                        'normal' : NORMAL['right']})
    if not block['front'] and not array.mask[px, y+1, z]:
        # Back pixel exist... create the back face.
        color = color_map.to_rgba(array[px, y+1, z])
        pixels.append({'vertex_coord' : [[pos[0], pos_u[1], pos[2]],
                                          [pos_u[0], pos_u[1], pos[2]],
                                          [pos[0], pos_u[1], pos_u[2]],
                                          [pos_u[0], pos_u[1], pos_u[2]]],
                             'color' : [color, color, color, color],
                        'normal' : NORMAL['back']
                                })
    if not block['left'] and not array.mask[px-1, y, z]:
        # Left pixel exist... create the left face.
        color = color_map.to_rgba(array[px - 1, y, z])
        pixels.append({'vertex_coord' : [[pos[0], pos[1], pos[2]],
                                         [pos[0], pos_u[1], pos[2]],
                                         [pos[0], pos[1], pos_u[2]],
                                         [pos[0], pos_u[1], pos_u[2]]],
                        'color' : [color, color, color, color],
                        'normal' : NORMAL['left']})
    if not block['back'] and not array.mask[px, y-1, z]:
        # Front pixel exist... create the fron face.
        color = color_map.to_rgba(array[px, y-1, z])
        pixels.append({'vertex_coord' : [[pos[0], pos[1], pos[2]],
                                          [pos_u[0], pos[1], pos[2]],
                                          [pos[0], pos[1], pos_u[2]],
                                          [pos_u[0], pos[1], pos_u[2]]],
                       'color' : [color, color, color, color],
                        'normal' : NORMAL['front']
                                })
    if not block['top'] and not array.mask[px, y, z+1]:
        # Top pixel exist... create the top face.
        color = color_map.to_rgba(array[px, y, z+1])
        pixels.append({'vertex_coord' : [[pos[0], pos[1], pos_u[2]],
                                              [pos[0], pos_u[1], pos_u[2]],
                                              [pos_u[0], pos[1], pos_u[2]],
                                              [pos_u[0], pos_u[1], pos_u[2]]],
                            'color' : [color, color, color, color],
                        'normal' : NORMAL['top']})
    return pixels
def generate_vertex(array, px, y, z, color_map, face = 'floor'):
    """
        Generate the vertex in the defined position.
    """
    voxel = None
    if not ismasked(array[px, y, z]):
#    if type(array[px, y, z]) not in [numpy.ma.core.MaskedArray,  # old numpy
#                                     numpy.ma.core.MaskedConstant]:
        factor = float(max(array.shape))
        pos = [((float(px) / factor) - array.shape[0] / (2 * factor)),
               ((float(y) / factor) - array.shape[1] / (2 * factor)),
               ((float(z) / factor) - array.shape[2] / (2 * factor))]
        pos_u = [((float(px + 1) / factor) - array.shape[0] / (2 * factor)),
               ((float(y+1) / factor) - array.shape[1] / (2 * factor)),
               ((float(z+1) / factor) - array.shape[2] / (2 * factor))]
               
        # The vertex is a 3d vertex composed on 6 faces
        #
        #       (x, y+1, z+1)_______
        #                  /·     /| (x+1, y+1, z+1)
        #                 / ·    / |
        #    (x, y, z+1) +------+  |
        #                |  ·······+ (x+1, y+1, z)
        #   (x, y+1, z)  | ·    | /      
        #                |·     |/       
        #                +------+ (x+1, y, z)
        #                (x, y, z)
        # By convention, the face of the bottom has the color form the pixel x,
        # y, z and the face from the top has the color of x, y, z+1.
        # Then, if pixel in z+1 is masked, we do not render it.
        # In order to not make twice the same face, we only consiter faces at
        # bottom : [(x, y, z), (x+1, y, z), (x+1, y+1, z), (x, y+1, z)]
        # front : [(x, y, z), (x+1, y, z), (x, y, z+1), (x+1, y, z+1)]
        # right : [(x, y, z), (x, y+1, z), (x, y, z+1), (x, y+1, z+1),]
        color = color_map.to_rgba(array[px, y, z])
        if face == 'floor' :
            if z+1 == array.shape[2]:
                pos[2] = pos_u[2]
            voxel = {'vertex_coord' : [[pos[0], pos[1], pos[2]],
                                              [pos[0], pos_u[1], pos[2]],
                                              [pos_u[0], pos[1], pos[2]],
                                              [pos_u[0], pos_u[1], pos[2]]],
                            'color' : [color, color, color, color],
                        'normal' : NORMAL['bottom'],
                            }
        elif face == 'front':
            if y+1 == array.shape[1]:
                pos[1] = pos_u[1]
            voxel = {'vertex_coord' : [[pos[0], pos[1], pos[2]],
                                          [pos_u[0], pos[1], pos[2]],
                                          [pos[0], pos[1], pos_u[2]],
                                          [pos_u[0], pos[1], pos_u[2]]],
                            'color' : [color, color, color, color],
                        'normal' : NORMAL['front']
                            }
        elif face == 'right':
            if px + 1 == array.shape[0]:
                pos[0] = pos_u[0]
            voxel = {'vertex_coord' : [[pos[0], pos[1], pos[2]],
                                          [pos[0], pos_u[1], pos[2]],
                                          [pos[0], pos[1], pos_u[2]],
                                          [pos[0], pos_u[1], pos_u[2]]],
                            'color' : [color, color, color, color],
                        'normal' : NORMAL['right']}
    return voxel
def make_coordinate(array, color_map):
    """
        This function creates the coordinate of the vertex and the corresponding
        color.
        
        It returns a list of dictionnary.
    """
    #array = array / array.max() * 4
#    color_array = color_array / color_array.max()
#    plot_array = pcolor(color_array)
    coordinate = list()
    factor = float(max(array.shape))
    for x in range(array.shape[0]):
        for y in range(array.shape[1]):
            z = array[px, y] / 10.
            if type(color_array[px, y]) == numpy.ma.core.MaskedArray:
                color = [1,1,1]
            else:
                color = color_map.to_rgba(color_array[px,y])
            pos_x = float(x) / array.shape[0]
            pos_y = float(y) / array.shape[1]
            pos_x_u = float(x + 1) / array.shape[0]
            pos_y_u = float(y + 1) / array.shape[1]
            this_pixel = {'tex_coord' : [[pos_x, pos_y],
                                      [pos_x, pos_y_u],
                                      [pos_x_u, pos_y],
                                      [pos_x_u, pos_y_u]],
                          'vertex_coord' : [[px / factor - 0.5,
                                             y / factor - 0.5,
                                             z], # bottom left
                                            [px / factor - 0.5,
                                             (y + 1)/factor - 0.5,
                                             z], # upper left
                                            [(px + 1)/factor - 0.5,
                                             y/factor - 0.5,
                                             z], # bottom right
                                            [(px + 1)/factor - 0.5,
                                             (y + 1)/factor - 0.5,
                                             z]], # upper right
                          'color' : [color, color, color, color]}
            coordinate.append(this_pixel)
            if x > 0:
                z_prev = array[px,-1, y] / 10
                # we make a step pixel in x.
                # The x and y are like pos_x and pos_y, pos_y_u
                # The pos_z are from this pixel and from the next.
                if type(color_array[px,-1,y]) == numpy.ma.core.MaskedArray:
                    color_prev = [1,1,1]
                else:
                    color_prev = color_map.to_rgba(color_array[px,-1,y])
                this_step = {'tex_coord' : [[pos_x, pos_y],
                                      [pos_x, pos_y_u],
                                      [pos_x, pos_y],
                                      [pos_x, pos_y_u]],
                          'vertex_coord' : [[px / factor - 0.5,
                                             y / factor - 0.5, 
                                             z], # bottom left
                                            [px / factor - 0.5,
                                             (y + 1) / factor - 0.5, 
                                             z], # upper left
                                            [px / factor - 0.5,
                                             y/factor - 0.5, 
                                             z_prev], # bottom right
                                            [px / factor - 0.5,
                                             (y + 1) / factor - 0.5, 
                                             z_prev]], # upper right
                          'color' : [color, color, color_prev, color_prev]}
                coordinate.append(this_step)
            if y > 0:
                z_prev = array[px, y-1] / 10
                # we make a step pixel in x.
                # The x and y are like pos_x and pos_y, pos_y_u
                # The pos_z are from this pixel and from the next.
                if type(color_array[px,y-1]) == numpy.ma.core.MaskedArray:
                    color_prev = [1,1,1]
                else:
                    color_prev = color_map.to_rgba(color_array[px,y-1])
                this_step = {'tex_coord' : [[pos_x, pos_y],
                                      [pos_x, pos_y],
                                      [pos_x_u, pos_y],
                                      [pos_x_u, pos_y]],
                          'vertex_coord' : [[px / factor - 0.5,
                                             y / factor - 0.5, 
                                             z], # bottom left
                                            [px / factor - 0.5,
                                             y / factor - 0.5, 
                                             z_prev], # upper left
                                            [(px + 1)/factor - 0.5,
                                             y / factor - 0.5, 
                                             z], # bottom right
                                            [(px + 1)/factor - 0.5,
                                             y / factor - 0.5, 
                                             z_prev]], # upper right
                          'color' : [color, color_prev, color, color_prev]}
                coordinate.append(this_step)
                                      
    return coordinate
NORMAL = {
    'top' : (0, 0, -1),
    'bottom' : (0, 0, 1),
    'left' : (0, -1, 0),
    'right' : (0, 1, 0),
    'front' : (-1, 0, 0),
    'back' : (1, 0, 0)
}
# Some api in the chain is translating the keystrokes to this octal string
# so instead of saying: ESCAPE = 27, we use the following.
ESCAPE = '\033'

if __name__ == '__main__':
    import numpy
    SIZE = 32
    ARRAY = numpy.random.randn(SIZE, SIZE)
    TMP_ARRAY = numpy.random.randn(SIZE, SIZE)
    COLOR_ARRAY = numpy.ma.array(TMP_ARRAY, mask = TMP_ARRAY < 0)
    TEST = Topography(array, COLOR_ARRAY)
