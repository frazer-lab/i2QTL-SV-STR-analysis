
# coding: utf-8

# In[2]:

from __future__ import division

import numpy as np
import os
import sys
import datetime
from subprocess import call
import subprocess
import glob
import djPyBio as DJ
import pandas as pd
import csv
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy 
import pybedtools
pd.set_option('display.max_columns', 500)


# In[13]:

# def format_labels_test(xlabel='none', ylabel='none', fontsizes=14, ticklabel_size = 14, 
#                   fontsize_xticklabel = False, fontsize_yticklabel = False, fontsize_xlabel = False, 
#                   fontsize_ylabel = False, rotation_x = 0, rotation_y=0):
    
#     print locals()


# In[14]:

# format_labels_test()


# In[33]:

def get_legend_patches(labels, just_rect=False, just_lines= False, linecolor= 'black', linewidth=0, marker ='o', 
                           markeredgecolor='black', markeredgewidth=1, 
                           markerfacecolor = 'white', markersize=20, linestyle= '-',
                       rectedgecolor='black', rectfacecolor = 'white', rectlinewidth=1, hatch = '', alpha = 1):
                        
    """Function to return a list of patches with any combination of markeredgewidths, facecolors, markersizes, edgecolors- only for points and other markers- lines with or without points
    
    if lines = True- return patches of just lines from colors_lines, and labels_lines
    to insert spaces in a legend- pass a list with linewidths or markersize= 0 and label = '' for those elements
    make a legend by ax.legend(handles =patches) and adjust as desired
    
    """

    patches = []
    num_patches = len(labels)
  
    props = [markeredgecolor, markerfacecolor, markeredgewidth, markersize, marker, 
             linewidth, linecolor, linestyle, rectedgecolor, rectfacecolor,  rectlinewidth, hatch, alpha]
    prop_names = ['markeredgecolor', 'markerfacecolor', 'markeredgewidth', 'markersize', 'marker', 
                  'linewidth', 'linecolor', 'linestyle','rectedgecolor', 'rectfacecolor', 
                  'rectlinewidth', 'hatch', 'alpha']
    

    for i, z in enumerate(props):
        if hasattr(z, '__iter__'): 
            assert len(z) == num_patches, 'incorrect length of list of {}'.format(prop_names[i])

        else:
            # duplicate the default into a list of length num patches- simplify patch list building below
            z = [z for l in range(num_patches)]
            props[i] = z

    (markeredgecolor, markerfacecolor, markeredgewidth, markersize, 
     marker, linewidth, linecolor, linestyle, rectedgecolor, rectfacecolor,  
     rectlinewidth, hatch, alpha) = props
    
    
    if just_lines:
        for i, l in enumerate(labels):
            # plot a line with no markers
            p= plt.Line2D((0,0),(0,0), linewidth=linewidth[i], color=linecolor[i], label = l, linestyle=linestyle[i])

            patches.append(p)     
    
    elif just_rect:
        for i, l in enumerate(labels):
            # plot rects not lines
            p = plt.Rectangle((0, 0), 0, 0, fc= rectfacecolor[i], ec = rectedgecolor[i], 
                              lw = rectlinewidth[i], linestyle = linestyle[i], hatch=hatch[i], label = l, alpha = alpha[i])
            patches.append(p)
    else:    
        for i, l in enumerate(labels):
            # plot lines/markers- default to no line just marker
            p= plt.Line2D((0,0),(0,0), linewidth= linewidth[i], color= linecolor[i], 
                          markeredgecolor = markeredgecolor[i],
                          markeredgewidth= markeredgewidth[i], markerfacecolor = markerfacecolor[i], 
                          markersize=markersize[i], label= l, marker= marker[i])
            patches.append(p)
        
    return patches

class label_annotator(object):
    """ experimental class used to define and apply styles to plots, allowing slight modifications on a per plot basis
    alternative to making a bunch of custom rcparams you edit"""
    def __init__(self, fontsizes=14, ticklabel_size = 14, 
                  fontsize_xticklabel = False, fontsize_yticklabel = False, fontsize_xlabel = False, 
                  fontsize_ylabel = False, rotation_x = 0, rotation_y=0):
        self.fontsizes= fontsizes
        self.fontsize_xticklabel = fontsize_xticklabel
        self.fontsize_yticklabel = fontsize_yticklabel
        self.fontsize_xlabel = fontsize_xlabel
        self.fontsize_ylabel = fontsize_ylabel
        self.rotation_x = rotation_x
        self.rotation_y  = rotation_y
        
    def format_labels(self, ax,  xlabel='none', ylabel='none'):

        kw_tls = [self.fontsize_xticklabel, self.fontsize_yticklabel]
        kw_axlabs = [self.fontsize_xlabel, self.fontsize_ylabel]

        ind_axis_labs = False
        for l in kw_axlabs:
            if l <> False:
                ind_axis_labs =True


        ind_tlax_labs = False
        for l in kw_tls:
            if l <> False:
                ind_tlax_labs = True


        if ind_tlax_labs == True:
            for k, axis in zip(kw_tls, ['x', 'y']):
                if k <> False:
                    ax.tick_params(axis= axis, labelsize = k)
                else:
                    ax.tick_params(axis=axis, labelsize = self.ticklabel_size)
        else:
            ax.tick_params(axis='both', labelsize = self.ticklabel_size)




        if ind_axis_labs == True:

            if self.fontsize_xlabel <> False:
                if xlabel <> 'none':
                    ax.set_xlabel(xlabel, fontsize=self.fontsize_xlabel)
                else:
                    ax.xaxis.get_label().set_fontsize(self.fontsize_xlabel)

                if self.fontsize_ylabel == False:
                    if ylabel <> 'none':
                        ax.set_ylabel(ylabel, fontsize=self.fontsizes)
                    else:
                        ax.yaxis.get_label().set_fontsize(self.fontsizes)



            if self.fontsize_ylabel <> False:
                if ylabel <> 'none':
                    ax.set_ylabel(ylabel, fontsize=self.fontsize_ylabel)
                else:
                    ax.yaxis.get_label().set_fontsize(self.fontsize_ylabel)

                if self.fontsize_xlabel == False:
                    if xlabel <> 'none':
                        ax.set_xlabel(xlabel, fontsize=self.fontsizes)
                    else:
                        ax.xaxis.get_label().set_fontsize(self.fontsizes)
        else:

            if xlabel <> 'none':
                ax.set_xlabel(xlabel, fontsize=self.fontsizes)
            else:
                ax.xaxis.get_label().set_fontsize(self.fontsizes)

            if ylabel <> 'none':
                ax.set_ylabel(ylabel, fontsize=self.fontsizes)
            else:
                ax.yaxis.get_label().set_fontsize(self.fontsizes)



        x_ticks = ax.get_xticklabels() 
        y_ticks = ax.get_yticklabels()

        for t in x_ticks:

            t.set_rotation(self.rotation_x)
        for t in y_ticks:
            t.set_rotation(self.rotation_y)




# In[ ]:

#  def format_labels(label_annotator_test, ax, fontsizes = get_val(fontsizes), ticklabel_size = get_val(ticklabel_size),
#                       fontsize_xticklabel = get_val(fontsize_xticklabel), fontsize_yticklabel = get_val(fontsize_yticklabel), 
#                       fontsize_xlabel = get_val(fontsize_xlabel), fontsize_ylabel = get_val(fontsize_ylabel), 
#                       rotation_x = get_val(rotation_x), rotation_y=get_val(rotation_y), 
#                       xlabel= get_val(xlabel), ylabel= get_val(ylabel):


# In[38]:

2+2


# In[40]:

initial_defaults = dict(fontsizes=14, ticklabel_size = 14, 
                  fontsize_xticklabel = False, fontsize_yticklabel = False, fontsize_xlabel = False, 
                  fontsize_ylabel = False, rotation_x = 0, rotation_y=0, xlabel= 'none', ylabel = 'none')


# In[51]:

class label_annotator_test(object):
    
    initial_defaults =  dict(fontsizes=14, ticklabel_size = 14,
                             fontsize_xticklabel = False, fontsize_yticklabel = False, fontsize_xlabel = False, 
                             fontsize_ylabel = False, rotation_x = 0, rotation_y=0, xlabel= 'none', ylabel = 'none')
    

    def __init__(self, fontsizes=14, ticklabel_size = 25, 
                  fontsize_xticklabel = False, fontsize_yticklabel = False, fontsize_xlabel = False, 
                  fontsize_ylabel = False, rotation_x = 0, rotation_y=0):
        self.fontsizes= fontsizes
        self.fontsize_xticklabel = fontsize_xticklabel
        self.fontsize_yticklabel = fontsize_yticklabel
        self.fontsize_xlabel = fontsize_xlabel
        self.fontsize_ylabel = fontsize_ylabel
        self.rotation_x = rotation_x
        self.rotation_y  = rotation_y
        
        self.annotator_defaults = locals()
        self.axes_params = locals()
     
        
        
        
        

    
    def format_labels(self, ax, xlabel = 'none',ylabel='none', fontsizes= None, ticklabel_size = None, 
                  fontsize_xticklabel = None, fontsize_yticklabel = None, fontsize_xlabel =None, 
                  fontsize_ylabel = None, rotation_x = None, rotation_y=None):
        
        d = locals()
        print d 
        for k, l in d.iteritems():
       
            if l <> None:
                self.axes_params[k] = l
                
        for k, l in self.axes_params.iteritems():
            setattr(self, k, l)
        
        self.change_plot()
        
    def change_plot(self):
        kw_tls = [self.fontsize_xticklabel, self.fontsize_yticklabel]
        kw_axlabs = [self.fontsize_xlabel, self.fontsize_ylabel]

        ind_axis_labs = False
        for l in kw_axlabs:
            if l <> False:
                ind_axis_labs =True
        
        ax = self.ax


        ind_tlax_labs = False
        for l in kw_tls:
            if l <> False:
                ind_tlax_labs = True


        if ind_tlax_labs == True:
            for k, axis in zip(kw_tls, ['x', 'y']):
                if k <> False:
                    ax.tick_params(axis= axis, labelsize = k)
                else:
                    ax.tick_params(axis=axis, labelsize = self.ticklabel_size)
        else:
            ax.tick_params(axis='both', labelsize = self.ticklabel_size)




        if ind_axis_labs == True:

            if self.fontsize_xlabel <> False:
                if self.xlabel <> 'none':
                    ax.set_xlabel(self.xlabel, fontsize=self.fontsize_xlabel)
                else:
                    ax.xaxis.get_label().set_fontsize(self.fontsize_xlabel)

                if self.fontsize_ylabel == False:
                    if self.ylabel <> 'none':
                        ax.set_ylabel(self.ylabel, fontsize=self.fontsizes)
                    else:
                        ax.yaxis.get_label().set_fontsize(self.fontsizes)



            if self.fontsize_ylabel <> False:
                if ylabel <> 'none':
                    ax.set_ylabel(ylabel, fontsize=self.fontsize_ylabel)
                else:
                    ax.yaxis.get_label().set_fontsize(self.fontsize_ylabel)

                if self.fontsize_xlabel == False:
                    if xlabel <> 'none':
                        ax.set_xlabel(self.xlabel, fontsize=self.fontsizes)
                    else:
                        ax.xaxis.get_label().set_fontsize(self.fontsizes)
        else:

            if self.xlabel <> 'none':
                ax.set_xlabel(self.xlabel, fontsize=self.fontsizes)
            else:
                ax.xaxis.get_label().set_fontsize(self.fontsizes)

            if self.ylabel <> 'none':
                ax.set_ylabel(self.ylabel, fontsize=self.fontsizes)
            else:
                ax.yaxis.get_label().set_fontsize(self.fontsizes)



        x_ticks = ax.get_xticklabels() 
        y_ticks = ax.get_yticklabels()

        for t in x_ticks:

            t.set_rotation(self.rotation_x)
        for t in y_ticks:
            t.set_rotation(self.rotation_y)
            
    def return_ax_params(self):
        print self.axes_params
            
    def return_annotator_defaults(self):
        print self.annotator_defaults
        





# In[ ]:




# In[15]:

def format_labels(ax, xlabel='none', ylabel='none', fontsizes=14, ticklabel_size = 14, 
                  fontsize_xticklabel = False, fontsize_yticklabel = False, fontsize_xlabel = False, 
                  fontsize_ylabel = False, rotation_x = 0, rotation_y=0):

    
    kw_tls = [fontsize_xticklabel, fontsize_yticklabel]
    kw_axlabs = [fontsize_xlabel, fontsize_ylabel]
    
    ind_axis_labs = False
    for l in kw_axlabs:
        if l <> False:
            ind_axis_labs =True
    
    
    ind_tlax_labs = False
    for l in kw_tls:
        if l <> False:
            ind_tlax_labs = True
            
            
    if ind_tlax_labs == True:
        for k, axis in zip(kw_tls, ['x', 'y']):
            if k <> False:
                ax.tick_params(axis= axis, labelsize = k)
            else:
                ax.tick_params(axis=axis, labelsize = ticklabel_size)
    else:
        ax.tick_params(axis='both', labelsize = ticklabel_size)
    
    
    
    
    if ind_axis_labs == True:
    
        if fontsize_xlabel <> False:
            if xlabel <> 'none':
                ax.set_xlabel(xlabel, fontsize=fontsize_xlabel)
            else:
                ax.xaxis.get_label().set_fontsize(fontsize_xlabel)

            if fontsize_ylabel == False:
                if ylabel <> 'none':
                    ax.set_ylabel(ylabel, fontsize=fontsizes)
                else:
                    ax.yaxis.get_label().set_fontsize(fontsizes)



        if fontsize_ylabel <> False:
            if ylabel <> 'none':
                ax.set_ylabel(ylabel, fontsize=fontsize_ylabel)
            else:
                ax.yaxis.get_label().set_fontsize(fontsize_ylabel)

            if fontsize_xlabel == False:
                if xlabel <> 'none':
                    ax.set_xlabel(xlabel, fontsize=fontsizes)
                else:
                    ax.xaxis.get_label().set_fontsize(fontsizes)
    else:
        
        if xlabel <> 'none':
            ax.set_xlabel(xlabel, fontsize=fontsizes)
        else:
            ax.xaxis.get_label().set_fontsize(fontsizes)

        if ylabel <> 'none':
            ax.set_ylabel(ylabel, fontsize=fontsizes)
        else:
            ax.yaxis.get_label().set_fontsize(fontsizes)
        
        
    
    x_ticks = ax.get_xticklabels() 
    y_ticks = ax.get_yticklabels()
    
    for t in x_ticks:
        
        t.set_rotation(rotation_x)
    for t in y_ticks:
        t.set_rotation(rotation_y)
        
    
    
        
        
        

def format_axes(ax, linewidth_gridlines=1, linewidth_spines = 3, color_gridlines = 'grey', color_spines = 'black', 
                     all_four=False, make_invisible=False):
    ticklines = ax.get_xticklines() + ax.get_yticklines()
    gridlines = ax.get_xgridlines() + ax.get_ygridlines()
    ticklabels = ax.get_xticklabels() + ax.get_yticklabels()

    for line in gridlines:
        line.set_color(color_gridlines)
        line.set_linewidth(linewidth_gridlines)
    
    for loc in ['bottom', 'left']:
        ax.spines[loc].set_color(color_spines)
        ax.spines[loc].set_linewidth(linewidth_spines)
    
    if all_four==True:
        for loc in ['right', 'top']:
            ax.spines[loc].set_color(color_spines)
            ax.spines[loc].set_linewidth(linewidth_spines)
            
            
    if make_invisible == True:
        ax.set_frame_on(False)
        ax.set_xlabel('')
        ax.set_xticklabels('')
        ax.set_yticklabels('')
        ax.grid(False)


# In[4]:

def add_text_to_ax(ax, text, loc = (-.1,1.1), font=False, fontsize = 18, **kwargs):
    
    if type(ax)== list:
        for a, lab in zip(ax, text):
            a.text(loc[0], loc[1], lab , ha='center', va='center', transform=a.transAxes, fontsize=fontsize, **kwargs)
    else:
        t = ax.text(loc[0], loc[1], text , ha='center', va='center', transform=ax.transAxes, fontsize=fontsize, **kwargs)
        return t


# In[5]:

def format_legend(ax, ncol=3, fontsize_legend = 14, title= '', fontsize_title = 18, bbox_to_anchor= (0,-.3, 0.8, 1),
                 fc = 'white', ec = 'black', ew = 3, frame_visible=False, legend_invisible=False, **kwargs):
    
    leg = ax.legend(title = title, bbox_to_anchor=bbox_to_anchor, fontsize = fontsize_legend, 
              loc=3,ncol=ncol , borderaxespad=0., frameon=True, **kwargs)
    
    if frame_visible == True:

        leg.get_frame().set_facecolor(fc)
        leg.get_frame().set_edgecolor(ec)
        leg.get_frame().set_linewidth(ew)
    if legend_invisible == True:
        leg.set_visible(False)
        
    if title <> '':
        ax.get_legend().get_title().set_fontsize(fontsize_title)
        
    
    return leg


def set_position(ax, left = 0, right = 0, bottom =0, top = 0):
    """ Set position of axes boundaries with exact coordinates (figure level) rather than width, height)"""
    
    width = right - left
    height = top - bottom
    
    ax.set_position([left, bottom, width, height])
    
  
def dict_gb_ranges(df, gb_cat):
    """ return the min and max coord of the elements of the bar plot after grouping them
    these define the boundaries of the bars below the plot
    
    df: dataframe pre-sorted in order of groupby [group1, group2]- (groups are columns)
    gb_cat: eg [group1, group2], how the data is sorted 
    items are expected to be plotted on integers from 0,l, 2- probably add a gap size option for plotting"""
    
    df['num'] = range(0,df.shape[0])
    dict_out = {}
    for i, d in df.groupby(gb_cat):
        min_out = d.num.min() -0.5
        max_out = d.num.max()+ 0.5
        dict_out[i] = [min_out, max_out]
    return dict_out


def annotate_region_below_axis(ax, left, right, bottom, height, fc = '#efe4bf',
                              ec = 'black', lw = 1, text='none', fontsize=15, zorder=14, boxstyle = 'round, pad=0',**kwargs):
    
    """ annotate boxes below an axis and annotate them with text- boxes are FancyBBox patch but could 
    also be drawn with brokenbarH objects, I don't know if there would be be better efficiency...
    but FancyBBox allows you to adjust the boxstyle etc (round, square) for better asthetic
    
    ax: ax to annotate
    left: left boundary desired (data units)
    right: right boundary desired (data units)
    bottom: bottom y coord of box (data units)
    height: height of box (in data units)"""
    
    
    width = right - left
    right = left + width
    top = bottom + height
    
    p = mpl.patches.FancyBboxPatch(
       (left, bottom), width, height,
        boxstyle=boxstyle, clip_on=False,
        fc= fc,
        ec = ec,  
        lw = lw, mutation_aspect = 1, zorder=zorder, **kwargs)
    ax.add_patch(p)
    
    ax.text(0.5*(left+right), 0.5*(bottom+top), text,
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=fontsize, color='black',
        transform=ax.transData, zorder=zorder+1)

def draw_bar_groupings(ax, d_ranges, color_dict, bottom, height, level = 0, label_on = True, fontsize = 10):
    """draw a bar for each grouping in a groupby that could be a single level or more
    ax: ax to annotate with bar
    d_ranges: the coords (x0,x1) of bar dimension
    color_dict: colors for the bars
    bottom: bottom of bars
    height: height of bars
    label_on: label the bars (True) with group name, or not
    fontsize: size of font on bars"""
    
    for l in d_ranges.keys():
        left, right = d_ranges[l]
        if level < 1:
            c = color_dict[l]
            label = l
        else:
            c = color_dict[l[level]]
            label = l[level]
        if not label_on:
            label = ''
            
        annotate_region_below_axis(ax, left, right, bottom, height, text=label, fc=c, fontsize=fontsize, zorder=13)

def annotate_multi_level(ax, df, group_order, color_dicts, bottoms, heights, label_bool, fontsizes=False):
    """ax: is the axes to annotate
    df: the df that will be used to plot the bars from.  This df has been presorted in the order of the groupby_order
    this means df = df.sort_values([group1, group2, group3]) 
    group_order: eg [group1, group2, group3]- sort order of the dataframe.  
    color_dicts: list of dictionaries of label:color for each group- in order of the group_order
    bottoms: list of bottom (y) level for the bars for each group- make this fit your ax as desired 
    heights: heights of each group- (top of each bar for each group will be at bottom+height)
    label_bool:  True/False for each group in group_order where True means label the boxes and false means don't label
    fontsizes: list of fontsizes, if not present- default 10 for all bars  """ 
    
    count = 1
    for g, cd, b, h, l in zip(group_order, color_dicts, bottoms, heights, label_bool):
        d_ranges = dict_gb_ranges(df, group_order[:count])
        if fontsizes: 
            fontsize = fontsizes[count-1]
            draw_bar_groupings(ax, d_ranges, cd, b, h, level = count - 1, label_on = l, fontsize=fontsize)
        else:
            draw_bar_groupings(ax, d_ranges, cd, b, h, level = count - 1, label_on = l)
            
        count +=1


def plot_per_sample_barplot(df, order_of_cols, colors, ax):
    


    df[order_of_cols].plot(kind='bar', stacked =True, ax = ax, colors = colors)



    ax.set_xticklabels('')
    ax.set_ylabel('number of SV sites')
    ax.tick_params(axis = 'x', top = False, bottom = False)
    ax.tick_params(axis = 'y', width = 1, length = 3)
    ax.xaxis.grid(False)

    leg = ax.legend(bbox_to_anchor = (1,1), loc = 'upper left', frameon=True)
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_edgecolor('black')

    group_order = ['STUDY', 'CELL_TYPE', 'ethnicity_mod']
    study_color_dict = {'HipSci': '#96b9f2', 'iPSCORE': '#95f1b0', 'HLI':'black'}
    ethnicity_color_dict = {'African': '#baced8', 'Non African': '#6b818c'}
    cell_type_color_dict = {'Blood': '#efe4bf', 'Fibroblast': '#dce8b2', 'iPSC':'#c2d6c8'}
    color_dict_order = [study_color_dict, cell_type_color_dict, ethnicity_color_dict]
    bottoms = [-200, -400, -600]
    heights = [200, 200, 200]
    label_bools = [True, True, False]


#     d_ranges_study = dict_gb_ranges(df, 'STUDY')

    divider = make_axes_locatable(ax)

    axsubx = divider.append_axes("bottom", size=0.8, pad=0.1)
    x_min, x_max = ax.get_xlim()
    axsubx.set_xlim(x_min, x_max)
    axsubx.set_ylim(-700,0)
    annotate_multi_level(axsubx, df, group_order, color_dict_order, bottoms, heights, label_bools)
    axsubx.set_axis_off()
    return ax





def annotate_with_boxes_and_text(ax, left, right, bottom, height, fc = '#efe4bf',
                                 ec = 'black', lw = 3, text='none', fontsize=15, zorder=14, **kwargs):
    
    """ place boxes with optional text on an axes given DATA coordinates"""
    width = right - left

    right = left + width
    top = bottom + height

    p = patches.FancyBboxPatch(
       (left, bottom), width, height,
        boxstyle='round, pad=0', clip_on=False,
        fc= fc,
        ec = ec,  
        lw = lw, mutation_aspect = 1, zorder=zorder, **kwargs)
    ax.add_patch(p)


    ax.text(0.5*(left+right), 0.5*(bottom+top), text,
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=fontsize, color='black',
        transform=ax.transData, zorder=zorder+1)

def fix_log_to_scalar_formatting(ax, axis = 'y'):
    """ take an axes in default log formatting from mpl- make the values scaler and remove the minor ticks"""
    
    formatter = mpl.ticker.ScalarFormatter()
    
    if axis in ['y', 'both']:
        
        ax.get_yaxis().set_major_formatter(formatter)
        ax.yaxis.get_major_formatter().set_scientific(False)
        ax.yaxis.get_minor_ticks()
        #             a.yaxis.get_major_formatter().set_useOffset(False)
        ax.get_yaxis().set_minor_formatter(formatter)
        ax.yaxis.get_minor_formatter().set_scientific(False)
    
    if axis in ['x', 'both']:
        ax.get_xaxis().set_major_formatter(formatter)
        ax.xaxis.get_major_formatter().set_scientific(False)
        ax.xaxis.get_minor_ticks()
        ax.get_xaxis().set_minor_formatter(formatter)
        ax.xaxis.get_minor_formatter().set_scientific(False)
    
    plt.sca(ax)
    plt.minorticks_off()
	
def format_base_pair_scale(x, log10 = False, mb_round = 2, kb_round = 1):
    """ input list of tick locations, output formatted ticks 
    with bp,kbpb, MB annotations - log: if the ticks
    are log scaled- what base is the log scale 
    example: log=10 for log base 10"""
    if log10:
        input_ticks = [(10**i) for i in x]
    else:
        input_ticks = x    
    out = []
    for t in input_ticks: 
        if abs(t) < 1000:
            f = "{}bp".format(t)
        
        elif (abs(t) >= 100000):
            f = "{}MB".format(round((t/1000000), mb_round))
            
        elif (abs(t) >= 1000):
            f = "{}kb".format(round((t/1000), kb_round))
            
        else:
            return 'missed case'
        out.append(f)
    
    return out

def label_offset_axes(ax, fig, text, x=-2, y  = 6, units = 'points', fontsize = 15, weight = 'bold', 
                      fontcolor= 'black'):
    """ Add text label offset from top left corner of some axis, certain  number of points or inches (see units) away
    from the corner 
    
    Note there are extra theoretically unnecessary conversions into figure units here- might not need this
    but this version works consistently, axes units did not- think I solved why but leaving like this, works"""
    inv = ax.transData.inverted()
    
    inv_fig = fig.transFigure.inverted()

    # specify the number of points or inches or dots to transform by
    # places an axes transform here for the axes on which you want to pick an anchor point on
    offsetter = mpl.transforms.offset_copy(ax.transAxes, fig=fig, x= x, y= y, units = 'points')

    # offset an anchor point - this will return display coordinates 
    corner_offset = offsetter.transform([0,1])
    
    fig_coord = inv_fig.transform(corner_offset)
    fig.text(fig_coord[0], fig_coord[1], text, fontsize = fontsize, weight = weight, color = 'black')


def format_easy_read_log(x, log10 = False, rounding = 1, ints = False):
    """ input list of tick locations, output formatted ticks 
    with bp,kbpb, MB annotations - log: if the ticks
    are log scaled- what base is the log scale 
    example: log=10 for log base 10"""
    if log10:
        input_ticks = [(10**i) for i in x]
    else:
        input_ticks = x    
    out = []
    for t in input_ticks: 
        if t < 1000:
            if not ints:
                f = copy.deepcopy(t)
            else:
                f = int(t)
        
        elif t > 100000:
            if not ints:
                f = "{}M".format(round((t/1000000), rounding))
            else:
                f = "{}M".format(int(round((t/1000000), rounding)))
            
        elif t >= 1000:
            if not ints:
                f = "{}K".format(round((t/1000), rounding))
            else:
                 f = "{}K".format(int(round((t/1000), rounding)))
            
        else:
            return 'missed case'
        out.append(f)
    
    return out


def as_si(x, ndp):
    """ format text as x10^power rather than scientific 1Epower format"""
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

