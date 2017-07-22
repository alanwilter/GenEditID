#!bin/python3

from plotly import tools
import plotly.offline as pyoff

dataplot = [
  dict(
    name = 'trace1',
    legendgroup = 'a',
    x = [1,2,3,4],
    y = [1,1,1,1],
    mode = 'markers',
    marker = {'symbol':'circle', 'color':'orange'}),
  dict(
    name = 'trace2',
    legendgroup = 'b',
    x = [1,2,3,4],
    y = [2,2,2,2],
    mode = 'markers',
    marker = {'symbol':'circle', 'color': 'steelblue'}),
  dict(
    name = 'trace1',
    legendgroup = 'a',
    showlegend = False,
    x = [1,2,3,4],
    y = [3,3,3,3],
    mode = 'markers',
    marker = {'symbol':'square', 'color':'orange'}),
  dict(
    name = 'trace2',
    legendgroup = 'b',
    showlegend = False,
    x = [1,2,3,4],
    y = [4,4,4,4],
    mode = 'markers',
    marker = {'symbol':'square', 'color':'steelblue'})
  ]

pyoff.plot(dict(data = dataplot))

# -----------------
# same with subplots

trace1a = dict(
   name = 'trace1',
   legendgroup = 'a',
   x = [1,2,3,4],
   y = [1,1,1,1],
   mode = 'markers',
   marker = {'symbol':'circle', 'color':'orange'})
trace1b = dict(
   name = 'trace2',
   legendgroup = 'b',
   x = [1,2,3,4],
   y = [2,2,2,2],
   mode = 'markers',
   marker = {'symbol':'circle', 'color': 'steelblue'})

trace2a = dict(
   name = 'trace1',
   legendgroup = 'a',
   showlegend = False,
   x = [1,2,3,4],
   y = [3,3,3,3],
   mode = 'markers',
   marker = {'symbol':'square', 'color':'orange'})
trace2b = dict(
   name = 'trace2',
   legendgroup = 'b',
   showlegend = False,
   x = [1,2,3,4],
   y = [4,4,4,4],
   mode = 'markers',
   marker = {'symbol':'square', 'color':'steelblue'})


fig = tools.make_subplots(rows=1, cols=2)
fig.append_trace(trace1a, 1, 1)
fig.append_trace(trace1b, 1, 1)
fig.append_trace(trace2a, 1, 2)
fig.append_trace(trace2b, 1, 2)

pyoff.plot(fig)
   