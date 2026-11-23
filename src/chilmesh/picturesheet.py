"""Deterministic visual summaries of planar mesh checks."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, Polygon, Rectangle, Patch
from matplotlib.lines import Line2D

from .scoresheet import _input


def draw_sheet(mesh_or_arrays, sheet, out_path, *, title, subtitle_lines, zoom_windows=None, region=None, partial=False):
    if isinstance(mesh_or_arrays, tuple):
        points, elements = _input(*mesh_or_arrays)
    else:
        points, elements = _input(mesh_or_arrays)
    polygons = [points[e] for e in elements]
    qualities = np.array([m['edge_angle_quality'] for m in sheet.element_measures])
    ratios = np.array([m['size_ratio'] for m in sheet.element_measures])
    finite = np.flatnonzero(np.isfinite(qualities)); worst = finite[np.argsort(qualities[finite], kind='stable')]
    bounds = np.ptp(points, axis=0) if len(points) else np.ones(2)
    scale = max(float(np.max(bounds)), 1e-12)
    if zoom_windows is None:
        defect_points = [np.mean(polygons[i], axis=0) for i in sheet.defects['inverted']['ids']]
        defect_points += [np.mean(polygons[i], axis=0) for i in sheet.defects['nonconvex']['ids']]
        defect_points += [np.mean(np.vstack([polygons[i], polygons[j]]), axis=0) for i,j in sheet.defects['overlap']['ids']]
        defect_points += [points[i] for i in sheet.defects['hanging_node']['ids']]
        defect_points += [points[i] for pair in sheet.defects['duplicate_vertex']['ids'] for i in pair]
        defect_points += [np.mean(poly, axis=0) for poly in sheet.hole_polygons]
        radius = max(scale*.09, 1e-12)
        clusters = []
        for point in defect_points:
            for cluster in clusters:
                proposed = np.asarray(cluster + [point])
                if np.all(np.ptp(proposed, axis=0) <= 2*radius):
                    cluster.append(point)
                    break
            else:
                clusters.append([point])
        candidates = [(np.min(cluster,axis=0)+np.max(cluster,axis=0))/2 for cluster in clusters]
        candidates += [np.mean(polygons[i],axis=0) for i in worst]
        candidates += [np.mean(p,axis=0) for p in polygons]
        zoom_windows = []
        for center in candidates:
            box = (center[0]-radius,center[1]-radius,center[0]+radius,center[1]+radius)
            if any(box[0] < other[2] and box[2] > other[0] and box[1] < other[3] and box[3] > other[1] for other in zoom_windows): continue
            zoom_windows.append(box)
            if len(zoom_windows) == 4: break
        for center in candidates:
            if len(zoom_windows) == 4: break
            zoom_windows.append((center[0]-radius,center[1]-radius,center[0]+radius,center[1]+radius))
        while len(zoom_windows)<4: zoom_windows.append((0,0,1,1))
    if len(zoom_windows) != 4: raise ValueError('zoom_windows must contain four boxes')
    fig = plt.figure(figsize=(16,12), dpi=150, layout='constrained')
    gs = fig.add_gridspec(5,4,height_ratios=[.5, 3, 3, 1.6, .18])
    header = fig.add_subplot(gs[0,:]); header.axis('off')
    header.text(0, .95, title, fontsize=17, weight='bold', va='top', transform=header.transAxes)
    header.text(0, .55, '\n'.join(subtitle_lines), fontsize=9, va='top', transform=header.transAxes)
    stamps = list(sheet.checks)
    if partial: stamps.append(dict(name='not complete',passed=False,detail=f'{sheet.holes} holes'))
    for i,check in enumerate(stamps):
        header.text(.99, .95-i*.37, f"{check['name']}: {'PASS' if check['passed'] else 'FAIL'}  {check['detail']}", transform=header.transAxes, ha='right', va='top', fontsize=8, bbox=dict(facecolor='white',edgecolor='#007c50' if check['passed'] else '#c43b35',boxstyle='round,pad=.3'))
    axa = fig.add_subplot(gs[1,:2]); axb = fig.add_subplot(gs[1,2:]); axc = fig.add_subplot(gs[2,:2]); zoomgrid = gs[2,2:].subgridspec(2,2)
    zoomaxes = [fig.add_subplot(zoomgrid[i//2,i%2]) for i in range(4)]
    def setup(ax, name, adjustable='datalim'):
        ax.set_title(name, loc='left', fontsize=10, weight='bold'); ax.set_aspect('equal', adjustable=adjustable); ax.autoscale_view(); ax.tick_params(labelsize=7)
    def base(ax):
        for p in sheet.hole_polygons: ax.add_patch(Polygon(p,facecolor='none',edgecolor='#bb3333',hatch='///',linewidth=1,zorder=0))
        inverted = set(sheet.defects['inverted']['ids'])
        ax.add_collection(PolyCollection(polygons, facecolors=['#c43b35' if i in inverted else '#56b4e9' if len(e)==4 else '#e69f00' for i,e in enumerate(elements)], edgecolors='#36454f', linewidths=.35, zorder=1))
        for i in sheet.defects['nonconvex']['ids']:
            ax.add_patch(Polygon(polygons[i], facecolor='#c43b35', edgecolor='#5b2080', hatch='xx', linewidth=1.4, zorder=2))
        for pair in sheet.defects['overlap']['ids']:
            for i in pair: ax.add_patch(Polygon(polygons[i], fill=False, edgecolor='#bb0000', linewidth=1.4, zorder=3))
        hanging = sheet.defects['hanging_node']['ids']
        duplicates = sorted({i for pair in sheet.defects['duplicate_vertex']['ids'] for i in pair})
        if hanging: ax.scatter(points[hanging,0],points[hanging,1],marker='x',c='#bb0000',s=35,zorder=5)
        if duplicates: ax.scatter(points[duplicates,0],points[duplicates,1],marker='s',facecolors='none',edgecolors='#bb0000',s=40,zorder=5)
        ax.autoscale()
    base(axa); setup(axa, 'A  Mesh')
    for key, marker in [('degree_3','o'),('degree_5_plus','x')]:
        ids = sheet.irregular_vertices[key]
        if ids: axa.scatter(points[ids,0],points[ids,1], marker=marker, c='#cc79a7', s=np.clip(400/np.sqrt(max(len(points),1)),2,25), alpha=.65, label=key, zorder=4)
    handles, labels = axa.get_legend_handles_labels()
    if sheet.defects['inverted']['count']: handles.append(Patch(facecolor='#c43b35',label='Inverted'))
    if sheet.defects['nonconvex']['count']: handles.append(Patch(facecolor='#c43b35',edgecolor='#5b2080',hatch='xx',label='Nonconvex'))
    if sheet.defects['overlap']['count']: handles.append(Patch(facecolor='none',edgecolor='#bb0000',label='Overlap'))
    if sheet.defects['hanging_node']['count']: handles.append(Line2D([],[],color='#bb0000',marker='x',linestyle='',label='Hanging node'))
    if sheet.defects['duplicate_vertex']['count']: handles.append(Line2D([],[],color='#bb0000',marker='s',fillstyle='none',linestyle='',label='Duplicate vertex'))
    if sheet.holes: handles.append(Patch(facecolor='none',edgecolor='#bb3333',hatch='///',label='Gap'))
    if handles: axa.legend(handles=handles,fontsize=7,loc='upper right')
    for i,(x0,y0,x1,y1) in enumerate(zoom_windows):
        axa.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,edgecolor='#d55e00',linewidth=1)); axa.text(x0,y1,str(i+1),color='#d55e00',fontsize=9,weight='bold')
    collection = PolyCollection(polygons,array=np.ma.masked_invalid(qualities),cmap='viridis',clim=(0,1),edgecolors='#444444',linewidths=.25)
    axb.add_collection(collection); axb.autoscale(); setup(axb,'B  Edge-angle quality'); fig.colorbar(collection,ax=axb,shrink=.75)
    poor = [i for i in worst if qualities[i] < 1 - 1e-9][:10]
    axb.text(.02,.02,f"Circled: {len(poor)}",transform=axb.transAxes,fontsize=8)
    for i in poor:
        center = np.mean(polygons[i],axis=0); axb.add_patch(Circle(center,scale*.012,fill=False,color='#d55e00',linewidth=1))
    if np.isfinite(ratios).any():
        extent = max(float(np.max(np.abs(ratios[np.isfinite(ratios)]-1))),.01)
        coll = PolyCollection(polygons,array=np.ma.masked_invalid(ratios),cmap='coolwarm',norm=TwoSlopeNorm(vmin=1-extent,vcenter=1,vmax=1+extent),edgecolors='#444444',linewidths=.25)
        axc.add_collection(coll); axc.autoscale(); fig.colorbar(coll,ax=axc,shrink=.75)
    else:
        axc.text(.5,.5,'No target size supplied',transform=axc.transAxes,ha='center'); axc.set_axis_off()
    setup(axc,'C  Size ratio')
    for i,(ax,box) in enumerate(zip(zoomaxes,zoom_windows)):
        base(ax); ax.set_xlim(box[0],box[2]); ax.set_ylim(box[1],box[3]); setup(ax, f'D{i+1}  Zoom', adjustable='box')
    a = fig.add_subplot(gs[3,:2]); b = fig.add_subplot(gs[3,2:])
    angles = [angle for m in sheet.element_measures for angle in m['angles']]
    angle_limit = 360 if any(angle > 180 for angle in angles) else 180
    a.hist(angles,bins=np.linspace(0,angle_limit,angle_limit//10+1),color='#56b4e9');
    for mark in (45,90,135): a.axvline(mark,color='#d55e00',linestyle='--',linewidth=.9)
    a.set_title('E  Interior angles'); a.set_xlabel('Degrees')
    b.hist(qualities[np.isfinite(qualities)],bins=np.linspace(0,1,21),color='#009e73'); b.set_title('E  Element quality'); b.set_xlabel('Quality')
    footer = fig.add_subplot(gs[4,:]); footer.axis('off')
    footer.text(.5,.5, f"Elements {sheet.counts['elements']}    Triangles {sheet.counts['triangles']}    Quad share {sheet.counts['quad_share']:.1%}    Irregular vertices {sheet.irregular_vertices['count']}    Holes {sheet.holes}",ha='center',va='center',fontsize=10,transform=footer.transAxes)
    fig.savefig(out_path,dpi=150); plt.close(fig)
