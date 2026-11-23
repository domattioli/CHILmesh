import json
import numpy as np
import pytest
from chilmesh import score_mesh


def grid(n=10):
    points = np.array([(x,y) for y in range(n+1) for x in range(n+1)],float)
    elements = [[y*(n+1)+x,y*(n+1)+x+1,(y+1)*(n+1)+x+1,(y+1)*(n+1)+x] for y in range(n) for x in range(n)]
    return points,elements


def test_grid():
    p,e = grid(); s = score_mesh(p,e,region=np.array([[0,0],[10,0],[10,10],[0,10]]))
    assert s.counts['elements'] == 100
    assert s.irregular_vertices['count'] == s.holes == 0
    assert s.clean and s.meshed_share == pytest.approx(1)
    assert all(m['edge_angle_quality'] == pytest.approx(1) and m['scaled_jacobian'] == pytest.approx(1) for m in s.element_measures)
    json.dumps(s.to_dict(), allow_nan=False)


@pytest.mark.parametrize('points, expected', [
    ([[0,0],[1,0],[1,1],[0,1]], (1,1,1,0)),
    ([[0,0],[2,0],[2,1],[0,1]], (.795,1,.632,0)),
    ([[0,0],[1,0],[1.5,3**.5/2],[.5,3**.5/2]], (.639,.866,.816,0)),
    ([[0,0],[2,0],[1.5,1],[.5,1]], (.653,.894,.784,.5)),
])
def test_hand(points,expected):
    s = score_mesh(np.array(points),[[0,1,2,3]])
    assert [s.element_measures[0][k] for k in ('edge_angle_quality','scaled_jacobian','stretch','taper')] == pytest.approx(expected,abs=1e-3)


def with_defects(names):
    p,e = grid()
    if 'inverted' in names: e[0] = e[0][::-1]
    if 'gap' in names: e[44] = None
    if 'duplicate_vertex' in names:
        p = np.vstack((p, p[0])); e[1][0] = len(p)-1
        # Element 1 does not use vertex 0: duplicate a vertex it does use.
        p[-1] = p[1]
    if 'hanging_node' in names:
        p = np.vstack((p, [[2.5,1]])); h = len(p)-1
        # Split the upper cell into three triangles along the new vertex.
        old = e[12]; a,b,c,d = old
        e[12] = [a,h,d]; e.extend([[h,b,c],[h,c,d]])
    if 'overlap' in names:
        p = np.vstack((p, [[.25,.25],[.75,.25],[.75,.75],[.25,.75]]))
        e.append(list(range(len(p)-4,len(p))))
    return p, [element for element in e if element is not None]


def damaged():
    return with_defects({'inverted','overlap','hanging_node','duplicate_vertex','gap'})


@pytest.mark.parametrize('name,expected', [
    ('inverted',[0]), ('overlap',[[0,100]]), ('hanging_node',[121]),
    ('duplicate_vertex',[[1,121]]), ('gap',[0]),
])
def test_single_defect(name,expected):
    p,e = with_defects({name})
    s = score_mesh(p,e,region=[[0,0],[10,0],[10,10],[0,10]])
    assert s.defects[name]['ids'] == expected
    assert all(not defect['count'] for key,defect in s.defects.items() if key != name)
    if name == 'inverted': assert s.element_measures[0]['edge_angle_quality'] == 0


def test_defects():
    p,e = damaged(); s=score_mesh(p,e,region=[[0,0],[10,0],[10,10],[0,10]])
    assert s.defects['inverted']['ids'] == [0]
    assert s.defects['overlap']['ids'] == [[0,101]]
    assert s.defects['hanging_node']['ids'] == [122]
    assert s.defects['duplicate_vertex']['ids'] == [[1,121]]
    assert s.defects['gap']['ids'] == [0]
    assert not s.clean


def test_nonconvex_has_zero_quality():
    p=np.array([[0,0],[2,0],[.5,.5],[0,2]])
    assert score_mesh(p,[[0,1,2,3]]).element_measures[0]['edge_angle_quality'] == 0


@pytest.mark.slow
@pytest.mark.skipif(__import__('os').environ.get('CHILMESH_SLOW') != '1', reason='set CHILMESH_SLOW=1')
def test_large_grid_timing():
    import time
    p,e=grid(142)
    start=time.perf_counter()
    sheet=score_mesh(p,e)
    elapsed=time.perf_counter()-start
    print(f'20,164 quads: {elapsed:.2f}s')
    assert sheet.clean and elapsed < 30


def test_mixed_angles():
    from chilmesh import Mesh
    p=np.array([[0,0],[1,0],[1,1],[0,1],[2,0],[2,1]],float)
    e=[[0,1,2,3],[1,4,2],[4,5,2]]
    mesh=Mesh(points=p,connectivity=np.array([[0,1,2,3],[1,4,2,1],[4,5,2,4]]))
    s=score_mesh(mesh); expected=mesh.interior_angles()
    for i,m in enumerate(s.element_measures):
        assert m['angles'] == pytest.approx(expected[i,:len(m['angles'])],abs=1e-6)


def test_disconnected_coverage_and_target():
    p=np.array([[1,1],[2,1],[2,2],[1,2]],float)
    s=score_mesh(p,[[0,1,2,3]],region=[[0,0],[3,0],[3,3],[0,3]],target_size=lambda xy: 2.)
    assert s.holes == 1
    assert s.meshed_share == pytest.approx(1/9)
    assert s.statistics['size_ratio']['mean'] == pytest.approx(.5)
    s.add_check('custom', True, 'ok')
    assert s.checks[-1]['passed']


def test_reflex_angle_and_nonconvex_defect():
    p = np.array([[0,0],[2,0],[.5,.5],[0,2]])
    sheet = score_mesh(p, [[0,1,2,3]])
    assert sheet.element_measures[0]['angles'][2] == pytest.approx(233.130102354156)
    assert sheet.defects['nonconvex']['ids'] == [0]
    assert all(not value['count'] for name,value in sheet.defects.items() if name != 'nonconvex')
    assert not sheet.clean


def test_nearly_straight_angle():
    p = np.array([[0,0],[1,0],[2,5e-11],[0,1]], float)
    sheet = score_mesh(p, [[0,1,2,3]])
    assert sheet.element_measures[0]['angles'][1] == pytest.approx(180 - np.degrees(5e-11), abs=1e-12)
    assert sheet.element_measures[0]['scaled_jacobian'] == pytest.approx(5e-11, abs=1e-15)
