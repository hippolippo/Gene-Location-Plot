import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
import numpy as np
import math
from typing import List
import json

LEFT = 0
RIGHT = 1

class Triangle:
    '''
    An Object representing one of the triangles for a `Graph`
    '''

    def __init__(self, position, filled, direction, color):
        '''
        Creates a `Triangle` object, representing one triangle that will be in a `Graph`

        Parameters:
        - position: location in Mb, position of the "flat end" of the triangle
        - filled: whether the triangle should be filled or transparent. (typically used as filled=gene, unfilled=pseudogene)
        - direction: should be either `LEFT` or `RIGHT` (0,1). Describes the direction that the "pointed end" of the triangle will point. (typically used as RIGHT= + strand, LEFT= - strand)
        - color: the color of the triangle as a string. The color should match the key used in your `Output` object when calling `finalize`
        '''
        self.position = position
        self.filled = filled
        self.direction = direction
        self.height_offset = 0
        self.color = color


class Graph:
    '''
    Object representing the graph for one chromosome
    '''

    def __init__(self, triangles: List[Triangle], length, centromere, name=None):
        '''
        Creates a `Graph` object, representing one chromosome

        Parameters:
        - triangles: a list of `Triangle` objects
        - length: the length of the chromosome in Mb
        - centromere: the position of the centromere in Mb
        - name: the label for this chromosome. e.g. `"Chromosome 1"`
        '''
        self.triangles: List[Triangle] = triangles
        self.plength = centromere
        self.qlength = length-centromere
        self.length = length
        self.name = name


class Output:
    '''
    Comprised of multiple `Graph` objects which each represent one chromosome

    add a graph using `add_graph` (populates bottom up), add the scale and key using `finalize` and render the plot using `visualize`

    Visual settings must be changed before `add_graph` is used and are
    Scale, MarkerHeight, MarkerSlope, TinyGap, SmallGap, LargeGap, BarHeight, FontSize, LabelSize.
    Typically the default values work well and these do not need to be modified
    '''

    Scale = 10
    MarkerHeight = 3
    MarkerSlope = 1
    TinyGap = 1.5
    SmallGap = .1*MarkerHeight
    LargeGap = 3*MarkerHeight
    BarHeight = 1.2*MarkerHeight
    FontSize = 1.4*MarkerHeight
    LabelSize = 20

    def __init__(self):
        '''
        Creates an empty plot, takes no parameters
        '''
        self.fig, self.ax = plt.subplots()
        self.graph_offset = 0
        self.this_offset = 0
        self.graphs = []
    
    def finalize(self, keyinfo=None, keypos=None):
        '''
        Adds the scale marker and key (if provided) to the plot

        Parameters:
        - keyinfo should be provided in a dictionary with keys being colors and values being the label. e.g. `{"#b326ff": "GRs", "orange": "IRs"}`
        - keypos is the position of the key. This parameter only matters if keyinfo is not `None`. if keypos is set to `None`, the position of the key will be automatically determined.
        '''
        # Draw Line
        x = self.graphs[0].length
        line = np.array([[x-self.Scale*3, -self.BarHeight*2.5], [x-self.Scale*2, -self.BarHeight*2.5],
                         [x-self.Scale*2, -self.BarHeight*2.5-.3], [x-self.Scale*3, -self.BarHeight*2.5-.3]])
        self.ax.add_patch(Polygon(line, color="black"))
        text = TextPath((x-self.Scale*2.5, -self.BarHeight*4), f"{self.Scale} Mb", prop=FontProperties(size=self.FontSize, weight=100))
        # Center the text
        text = TextPath((x-self.Scale*2.5-text.get_extents().width/2, -self.BarHeight*4), f"{self.Scale} Mb", prop=FontProperties(size=self.FontSize, weight=100))
        plt.gca().add_patch(PathPatch(text, color="black",fill=True))

        # Draw Key
        if keyinfo is not None:
            if keypos is None:
                keypos = (self.graphs[0].length+self.Scale*2, self.Scale*3)
            kx, ky = keypos
            s = self.Scale * .6
            yoff = 0
            # Draw each item in the key
            for color, text in keyinfo.items():
                square = np.array([[kx,ky+yoff],[kx+s,ky+yoff],[kx+s,ky+yoff+s],[kx,ky+yoff+s]])
                self.ax.add_patch(Polygon(square, color=color, fill=True))
                textpath = TextPath((kx+s*1.3, ky+yoff+s*.1), text, prop=FontProperties(size=s, weight=100))
                plt.gca().add_patch(PathPatch(textpath, color="black", fill=True))
                yoff -= s * 1.7

    def visualize(self):
        '''
        Draws the plot in a matplotlib window
        - If you want a scale marker and key, call `finalize` first
        '''
        self.ax.axis('off')
        self.ax.autoscale()
        self.ax.set_aspect('equal')
        plt.show()

    def add_graph(self, graph: Graph):
        '''
        adds a `Graph` object to the plot
        - graphs are added from bottom up
        '''
        zones = [[] for _ in range(math.ceil(graph.length/self.MarkerHeight))]
        for triangle in graph.triangles:
            zone = triangle.position // self.MarkerHeight
            possible_conflicts = sum(zones[math.floor(max(0,zone-2)):math.ceil(zone+1)], start=[])
            had_conflict = True
            while had_conflict:
                had_conflict = False
                for possible_conflict in possible_conflicts:
                    clearance = self._required_clearance(triangle, possible_conflict)
                    if clearance != 0:
                        had_conflict = True
                    triangle.height_offset += clearance
            zones[round(zone)].append(triangle)
            self._draw_triangle(triangle)
        self._draw_bars(graph)
        if graph.name != None:
            text = TextPath((-self.MarkerHeight, self.graph_offset-self.BarHeight*2-self.FontSize), graph.name, prop=FontProperties(size=self.FontSize, weight=100))
            plt.gca().add_patch(PathPatch(text, color="black",fill=True))
        self.graph_offset = self.this_offset + self.LabelSize
        self.this_offset = self.graph_offset
        self.graphs.append(graph)

            
    def _required_clearance(self, placing: Triangle, reference: Triangle):
        if reference.direction == LEFT:
            if placing.direction == RIGHT:
                # if they are close enough together, stack it on top
                if placing.position - reference.position < self.SmallGap and reference.height_offset-self.MarkerHeight <= placing.height_offset and reference.height_offset+self.MarkerHeight >= placing.height_offset:
                    return reference.height_offset + self.MarkerHeight - placing.height_offset
                else:
                    return 0
            else:
                # Both Left
                if placing.position - reference.position > self.MarkerHeight / 2 / self.MarkerSlope + self.SmallGap:
                    # Far enough apart horizontally
                    return 0
                # Minimum height to be above reference
                clearance_height = reference.height_offset + self.MarkerHeight - self.MarkerSlope * (placing.position - reference.position)
                if placing.height_offset >= clearance_height:
                    # already high enough
                    return 0
                # Maximum height to be below reference
                below_height = reference.height_offset - self.LargeGap + self.MarkerSlope * (placing.position - reference.position) - self.MarkerHeight
                if placing.height_offset <= below_height:
                    # low enough
                    return 0
                # Enough to be at clearance height
                return clearance_height - placing.height_offset
        else:
            # Reference is pointing right
            
            # Check if far enough horizontally
            if placing.direction == LEFT:
                if placing.position - reference.position > self.MarkerHeight / self.MarkerSlope + self.SmallGap:
                    # Far enough apart horizontally
                    return 0
            else:
                if placing.position - reference.position > self.MarkerHeight / 2 / self.MarkerSlope + self.SmallGap:
                    # Far enough apart horizontally
                    return 0
                
            # Minimum height to be above reference
            clearance_height = reference.height_offset + self.MarkerHeight - self.MarkerSlope * (placing.position - reference.position)
            if placing.direction == LEFT and abs(placing.position - self.MarkerHeight / 2 / self.MarkerSlope - reference.position) < self.MarkerHeight / 2 / self.MarkerSlope * .3:
                # Facing opposite directions and overlapping by more than 30%
                clearance_height += self.TinyGap
            if placing.height_offset >= clearance_height:
                # already high enough
                return 0
            # Maximum height to be below reference
            below_height = reference.height_offset + self.MarkerSlope * (placing.position - reference.position) - self.LargeGap - self.MarkerHeight
            if placing.height_offset <= below_height:
                # low enough
                return 0
            # Enough to be at clearance height
            return clearance_height - placing.height_offset
        
    def _draw_bars(self, graph: Graph):
        '''
        Draws the bars below each chromosome
        '''
        graph.triangles.sort(key = lambda t: t.position)
        y1 = self.graph_offset - self.MarkerHeight/4 - self.BarHeight
        y2 = self.graph_offset - self.MarkerHeight/4
        def get_shape(x1, x2):
            divot = self.BarHeight/5
            return np.array([
                [x1+divot, y1],
                [x2-divot, y1],
                [x2, y1+self.BarHeight/3],
                [x2, y2-self.BarHeight/3],
                [x2-divot, y2],
                [x1+divot, y2],
                [x1, y2-self.BarHeight/3],
                [x1, y1+self.BarHeight/3]
            ])
        bar1 = get_shape(0,graph.plength)
        bar2 = get_shape(graph.plength, graph.length)

        self.ax.add_patch(Polygon(bar1, color="lightgrey", fill=True))
        self.ax.add_patch(Polygon(bar1, color="black", fill=False))

        self.ax.add_patch(Polygon(bar2, color="grey", fill=True))
        self.ax.add_patch(Polygon(bar2, color="black", fill=False))
    
    def _draw_triangle(self, triangle: Triangle):
        direction = (-1 if triangle.direction == LEFT else 1)
        shape = np.array([[triangle.position, triangle.height_offset+self.graph_offset],
                             [triangle.position, triangle.height_offset+self.MarkerHeight+self.graph_offset],
                             [triangle.position + direction * (self.MarkerHeight/(2*self.MarkerSlope)), triangle.height_offset + self.MarkerHeight/2+self.graph_offset]
                            ])
        if triangle.height_offset + self.MarkerHeight + self.graph_offset > self.this_offset:
            self.this_offset = self.graph_offset + triangle.height_offset + self.MarkerHeight

        self.ax.add_patch(Polygon(shape, color=triangle.color, fill=triangle.filled))

