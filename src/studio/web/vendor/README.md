# Vendored libraries

Kept in the package so that the studio works offline, as on a compute node
with no outbound network. Each file keeps its own license header.

| File | Library | License | Source |
| --- | --- | --- | --- |
| `plotly.min.js` | plotly.js 2.35.2, the official *cartesian* partial bundle (bar, scatter, heatmap, histogram, box, contour) | MIT | `https://cdn.plot.ly/plotly-cartesian-2.35.2.min.js` |
| `highlight.min.js` | highlight.js 11.9.0, core | BSD-3-Clause | `https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js` |
| `julia.min.js` | highlight.js 11.9.0, Julia grammar | BSD-3-Clause | `https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/languages/julia.min.js` |

Plotly.js is the version MFH Studio vendors, in the bundle that holds the 2-D
traces this studio draws. To refresh a file, download it again from its source
with the new version number, and record the version here.
