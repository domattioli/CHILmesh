import json

import numpy as np

from chilmesh import audit


CELL = np.array([[0., 0.], [1, 0], [1, 1], [0, 1]])


class Toy:
    def __init__(self, lenient=False):
        self.lenient = lenient

    def reset(self, seed):
        self.placed = []
        self.steps = 0

    def domain(self):
        return CELL

    def front(self):
        return CELL if not self.placed else np.empty((0, 2))

    def elements(self):
        return self.placed

    def sample_action(self, rng):
        return rng.integers(0, 1, size=2)

    def candidate(self, action):
        return CELL + action

    def apply(self, action):
        accepted = not self.placed or self.lenient
        if accepted:
            self.placed.append(self.candidate(action))
        self.steps += 1
        return accepted, self.steps == (2 if self.lenient else 1)


def test_honest():
    report = audit(Toy(), 50)
    assert report.passed
    assert report.counts['agree'] == 50
    json.dumps(report.to_dict())


def test_checker_kwargs():
    report = audit(Toy(), 1, checker_kwargs={'angle_margin_deg': 100})
    assert report.disagreements[0]['checker_reason'] == 'nonconvex'
    assert report.to_dict()['checker_kwargs'] == {'angle_margin_deg': 100}
    assert audit(Toy(), 1).to_dict()['checker_kwargs'] == {}


def test_lenient():
    # Keep a front for the second proposal so geometry can assess overlap.
    class Lenient(Toy):
        def front(self):
            return CELL
    report = audit(Lenient(True), 1, seed=17)
    assert report.counts['disagree'] == 1
    assert report.disagreements[0]['trial_seed'] == 17
    assert report.disagreements[0]['step'] == 1


def test_candidate_matching_cyclic_order_and_scaled_tolerance():
    class Rotated(Toy):
        def domain(self):
            return CELL * 100

        def front(self):
            return self.domain() if not self.placed else np.empty((0, 2))

        def candidate(self, action):
            return self.domain()

        def apply(self, action):
            self.placed.append(np.roll(self.candidate(action), 2, axis=0) + 5e-8)
            return True, True

    assert audit(Rotated(), 1).passed

    class Reversed(Rotated):
        def apply(self, action):
            self.placed.append(self.candidate(action)[::-1])
            return True, True

    assert 'candidate_not_added' in audit(Reversed(), 1).invariant_failures[0]['failed_keys']
