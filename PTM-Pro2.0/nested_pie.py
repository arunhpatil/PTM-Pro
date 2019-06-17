"""
=================
Nested pie charts
=================

The following examples show two ways to build a nested pie chart
in Matplotlib. Such charts are often referred to as donut charts.

"""

import matplotlib.pyplot as plt
import numpy as np
from pie_inc import *
import sys
pwd=sys.argv[1]
###############################################################################
# The most straightforward way to build a pie chart is to use the
# :meth:`pie method <matplotlib.axes.Axes.pie>`
#
# In this case, pie takes values corresponding to counts in a group.
# We'll first generate some fake data, corresponding to three groups.
# In the inner circle, we'll treat each number as belonging to its
# own group. In the outer circle, we'll plot them as members of their
# original 3 groups.
#
# The effect of the donut shape is achieved by setting a `width` to
# the pie's wedges through the `wedgeprops` argument.


fig, ax = plt.subplots()


size = 0.3
vals = np.array(val_list)

cmap = plt.get_cmap("tab20c")
outer_colors = cmap(np.arange(4)*4)
inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))

wedges1, texts1 =ax.pie(vals.sum(axis=1), radius=1-size, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'))

wedges, texts = ax.pie(vals.flatten(), radius=1, colors=inner_colors,
       wedgeprops=dict(width=size, edgecolor='w'))

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=5, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    if unique_sites[i] != "0":
        ax.annotate(unique_sites[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),horizontalalignment=horizontalalignment, **kw)

ax.legend(wedges1, ptms,
          title="PTMS",
          loc="center left",
          bbox_to_anchor=(-0.8,0.5, 1, 0.5))

ax.set(aspect="equal", title='PTM-Pro2.0')
plt.show()
# Pad the saved area by 10% in the x-direction and 20% in the y-direction
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig(pwd+'/PTM_summary.png', bbox_inches=extent.expanded(2.1, 2.2))
fig.savefig(pwd+'/PTM_summary.eps', bbox_inches=extent.expanded(2.1, 2.2))
