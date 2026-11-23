import matplotlib.image as mpimg
import pytest
from chilmesh import score_mesh, draw_sheet
from test_scoresheet import grid, damaged

@pytest.mark.parametrize('factory,partial', [(grid,False),(damaged,True)])
def test_picture(tmp_path,factory,partial):
    p,e=factory(); sheet=score_mesh(p,e,region=[[0,0],[10,0],[10,10],[0,10]])
    path=tmp_path/'sheet.png'
    draw_sheet((p,e),sheet,path,title='Mesh',subtitle_lines=['Check'],partial=partial)
    assert path.stat().st_size > 10000
    assert mpimg.imread(path).shape[:2] == (1800,2400)


def test_reflex_histogram_and_nonconvex_legend(tmp_path, monkeypatch):
    import numpy as np
    from matplotlib.axes import Axes
    bins_seen = []
    legends = []
    original_hist, original_legend = Axes.hist, Axes.legend
    def capture_hist(self, *args, **kwargs):
        bins_seen.append(kwargs.get('bins'))
        return original_hist(self, *args, **kwargs)
    def capture_legend(self, *args, **kwargs):
        legends.extend(h.get_label() for h in kwargs.get('handles', []))
        return original_legend(self, *args, **kwargs)
    monkeypatch.setattr(Axes, 'hist', capture_hist)
    monkeypatch.setattr(Axes, 'legend', capture_legend)
    points = np.array([[0,0],[2,0],[.5,.5],[0,2]])
    elements = [[0,1,2,3]]
    sheet = score_mesh(points,elements)
    draw_sheet((points,elements),sheet,tmp_path/'reflex.png',title='Reflex',subtitle_lines=[])
    assert bins_seen[0][-1] == 360
    assert 'Nonconvex' in legends
