"""
========
Barchart
========

A bar plot with errorbars and height labels on individual bars.
"""
import numpy as np
import matplotlib.pyplot as plt
from bar_inc import *
import sys
pwd=sys.argv[1]


ind = np.arange(len(unambiguous_peptides))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, unambiguous_peptides, width,
                color='SkyBlue', label='Unambiguous')
rects2 = ax.bar(ind + width/2, ambiguous_peptides, width, 
                color='IndianRed', label='Ambiguous')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of peptides')
ax.set_title('Summary of modified peptides')
ax.set_xticks(ind)
ax.set_xticklabels(ptms)
ax.legend()


def autolabel(rects, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                '{}'.format(height), ha=ha[xpos], va='bottom')


autolabel(rects1)
autolabel(rects2)

#ax.set(aspect="equal", title='PTM-Pro2.0')
#plt.show()
# Pad the saved area by 10% in the x-direction and 20% in the y-direction
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig(pwd+'/peptide_ambiguity.png')
fig.savefig(pwd+'/peptide_ambiguity.eps')

