#!/usr/bin/env python
"""Plot data from a SLURM HDF5 profile file generated with sh5util"""

from __future__ import division, print_function

try:
    import h5py
except ImportError:
    print("Do 'module load Python' to get a version with HDF5 support")
    raise

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt, colors, gridspec
import sys, time, datetime, numpy

to_rgba = colors.ColorConverter().to_rgba

def the_only(seq):
    seq = iter(seq)
    item = next(seq)
    sentinal = object()
    assert next(seq, sentinal) is sentinal, "More than one"
    return item

def pick_io_scale(max_io):
    if max_io > 200:
        io_scale = 1000
        io_unit = 'GB'
    else:
        io_scale = 1
        io_unit = 'MB'
    return (io_scale, io_unit)

STEPS = sys.argv[2:]
fn = sys.argv[1]
f = h5py.File(fn)
steps = f['Steps']

# list(f['Steps']['0']['Nodes']['compute-b1-065']['Tasks']['0'])[0]
# fields = f.values()[0].values()[0].values()[0].values()[0].values()[0].values()[0].dtype.fields
# dtype([('ElapsedTime', '<u8'), ('EpochTime', '<u8'), 
# ('CPUFrequency', '<u8'), ('CPUTime', '<f8'), ('CPUUtilization', '<f8'), 
# ('RSS', '<u8'), ('VMSize', '<u8'), ('Pages', '<u8'), ('ReadMB', '<f8'), ('WriteMB', '<f8')])
# (30, 1450264708, 2701000, 23.62, 78.73333333333333, 5290260, 5490448, 0, 517.7017669677734, 0.01776123046875)

n_tasks = n_nodes = 0
steps_to_show = set()

for (step_name, step) in steps.items():
    if STEPS and step_name not in STEPS:
        continue
    #elif step_name == 'batch':
    #    continue
    steps_to_show.add(step_name)
    for node in step['Nodes'].values():
        n_nodes += 1
        for task in node.get('Tasks', {}).values():
            n_tasks += 1
            
low_alpha = 1.0 / n_tasks
alpha = 0.5 + 0.5 ** n_tasks

n_steps = len(steps_to_show)
if n_steps > 1:
    by_step = True
    n_colors = n_steps
else:
    by_step = False
    n_colors = n_tasks
    
cmap = plt.get_cmap('hsv', n_colors+1)

if by_step:
    fig = plt.figure(figsize=(11.69-1,8.27-1))
    gs = gridspec.GridSpec(4, 2, width_ratios=[10,1], height_ratios=[1,1,1,1])
    cpu_subplot = plt.subplot(gs[0,0])
    mem_subplot = plt.subplot(gs[1,0], sharex=cpu_subplot)
    iio_subplot = plt.subplot(gs[2,0], sharex=cpu_subplot)
    cio_subplot = plt.subplot(gs[3,0], sharex=cpu_subplot)
    footer = plt.subplot(gs[:,1])
    #fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    for ax in [cpu_subplot, mem_subplot, iio_subplot]:
        plt.setp(ax.get_xticklabels(), visible=False)
    footer.axis('off')
    #(fig, (cpu_subplot, mem_subplot, iio_subplot, cio_subplot, footer)) = plt.subplots(5, 1, sharex=True, squeeze=True)
else:
    (fig, (cpu_subplot, mem_subplot, iio_subplot, cio_subplot)) = plt.subplots(4, 1, sharex=True, squeeze=True, figsize=(11.69-1,8.27-1))
    footer = None
    
fig.suptitle(fn.split('/')[-1])
fig.subplots_adjust(hspace=0.15)

# Written bytes not actually negative:
import matplotlib.ticker
class NotNegFormatter(matplotlib.ticker.ScalarFormatter):
    def __call__(self, x, pos=None):
        return super(type(self), self).__call__(abs(x), pos)
iio_subplot.yaxis.set_major_formatter(NotNegFormatter())
iio_subplot.axhline(y=0, c='black')


base_start_time = None
max_iio = max_cio = 0
iio_plots = []
cio_plots = []

i_colour = 0
legend_data = []
for (i_step, (step_name, step)) in enumerate(steps.items()):
    if step_name not in steps_to_show:
        continue
    for (i_node, node) in enumerate(step['Nodes'].values()):
        for task in node.get('Tasks', {}).values():
            if by_step:
                i_colour = i_step
                series_name = step_name
            else:
                i_colour += 1
                series_name = None
            c = to_rgba(cmap(i_colour), alpha)

            data = task #list(task)
            if len(data) < 2:
                continue

            start_time = data[0, 'EpochTime']
            if base_start_time is None:
                base_start_time = start_time
                
            x = (data[:, 'EpochTime'] - base_start_time) / 60
            # time points not guaranteed to be evenly spaced.
            dx = numpy.diff(x)
            filter = dx >= numpy.median(dx)/2
            #filter[x[:-1]>110] = False  # minutes
            
            dx = dx[filter]
            x = x[numpy.concatenate([filter, [False]])]
            filter = numpy.concatenate([[False], filter])
            
            y = numpy.array(data[:,'CPUUtilization'])[filter]
            #print(x.shape, y.shape)
            y[y<0] = 0 # would be better to mask
            y[y>2400] = 2400 # temporary hack for job 61072225 glitch
            #print(x.shape, y.shape)
            cpu_subplot.plot(x, y/100, color=c, label=series_name)
            
            mem_subplot.plot(x, data[:,'RSS'][filter] / 1e6, color=c)
            
            io_series = [('ReadMB', '-', 1), ('WriteMB', ':', -1)]
            io_data = numpy.array([data[:,series_name][filter] for (series_name, pattern, sign) in io_series])
            io_data[io_data < 0] = 0  # Why do negative numbers occur here?
            cio_data = numpy.add.accumulate(io_data, axis=-1)
            iio_data = io_data / (dx*60)  # to make rates per-second. 
            # Postpone I/O plotting until all I/O data has been seen so as to pick GB vs MB scale.
            max_iio = max(max_iio, iio_data.max())
            max_cio = max(max_cio, cio_data.max())
            iio_plots.extend(zip(io_series, [c, c], [x, x], [dx, dx], iio_data))
            cio_plots.extend(zip(io_series, [c, c], [x, x], [dx, dx], cio_data))

(iio_scale, iio_unit) = pick_io_scale(max_iio)
iio_subplot.set_ylabel(iio_unit + '/s')
for ((series_name, pattern, sign), colour, x, dx, y) in iio_plots:
    b = y!=0
    iio_subplot.bar(x[b], y[b]*sign/iio_scale, dx[b], color=colour, linewidth=0)

(cio_scale, cio_unit) = pick_io_scale(max_cio)
cio_subplot.set_ylabel(cio_unit)
for ((series_name, pattern, sign), colour, x, dx, y) in cio_plots:
    cio_subplot.plot(x, y/cio_scale, color=colour, linestyle=pattern)

cpu_subplot.set_ylabel("CPUs")
mem_subplot.set_ylabel("GB")

plt.autoscale(axis='x', tight=True)
plt.xlabel("Minutes")

if footer:
    handles = []
    labels = []
    for (h, l) in zip(*cpu_subplot.get_legend_handles_labels()):
        if l not in labels:
            handles.append(h)
            labels.append(l)
    footer.legend(handles, labels, mode='expand', loc='upper left')

plt.savefig(fn.split('/')[-1].split('.')[0]+'_profile.png')








