"""Compare a proposal source with independent single-loop geometry checks."""
from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np

from .front_validity import _area, check_front_invariants, check_quad_against_front


class ProposalSource(Protocol):
    """Interface for an advancing-front proposal source."""
    def reset(self, seed: int) -> None: ...
    def domain(self) -> np.ndarray: ...
    def front(self) -> np.ndarray: ...
    def elements(self) -> list[np.ndarray]: ...
    def sample_action(self, rng: np.random.Generator) -> Any: ...
    def candidate(self, action: Any) -> np.ndarray | None: ...
    def apply(self, action: Any) -> tuple[bool, bool]: ...


def _json(value):
    if isinstance(value, np.ndarray):
        return _json(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    return value


@dataclass
class AuditReport:
    """Counts and reproducible records from an audit."""
    counts: dict = field(default_factory=lambda: dict.fromkeys(
        ('trials', 'steps', 'accepted', 'rejected', 'refused', 'agree',
         'disagree', 'invariant_failures'), 0))
    disagreements: list = field(default_factory=list)
    invariant_failures: list = field(default_factory=list)
    closing: list = field(default_factory=list)
    checker_kwargs: dict = field(default_factory=dict)

    @property
    def passed(self):
        return not self.disagreements and not self.invariant_failures

    def to_dict(self):
        return _json(dict(counts=self.counts, disagreements=self.disagreements,
                          invariant_failures=self.invariant_failures,
                          closing=self.closing, checker_kwargs=self.checker_kwargs, passed=self.passed))


def audit(source: ProposalSource, trials: int, seed: int = 0, max_steps: int = 10_000,
          checker_kwargs: dict | None = None) -> AuditReport:
    """Audit trials with independent seeded action streams and CCW checks."""
    report = AuditReport(checker_kwargs=dict(checker_kwargs or {}))
    for trial in range(trials):
        trial_seed = seed + trial
        source.reset(trial_seed)
        rng = np.random.default_rng(trial_seed)
        domain = np.array(source.domain(), dtype=float, copy=True)
        reverse = _area(domain) < 0
        orient = lambda p: np.asarray(p, dtype=float)[::-1].copy() if reverse else np.asarray(p, dtype=float).copy()
        domain = orient(domain)
        tolerance = 1e-9 * min(np.linalg.norm(b - a) for a, b in
                               zip(domain, np.roll(domain, -1, axis=0)))
        report.counts['trials'] += 1
        for step in range(max_steps):
            front = orient(source.front())
            elements = [orient(e) for e in source.elements()]
            action = source.sample_action(rng)
            candidate = source.candidate(action)
            candidate = None if candidate is None else orient(candidate)
            checked = None if candidate is None else check_quad_against_front(
                front, candidate, elements, domain, **report.checker_kwargs)
            accepted, done = source.apply(action)
            report.counts['steps'] += 1
            report.counts['accepted' if accepted else 'rejected'] += 1
            record = dict(trial=trial, trial_seed=trial_seed, step=step)
            if checked is None:
                report.counts['refused'] += 1
            elif bool(accepted) == checked[0]:
                report.counts['agree'] += 1
            else:
                report.counts['disagree'] += 1
                report.disagreements.append(dict(record, action=_json(action), candidate=candidate.tolist(),
                                                 front=front.tolist(), checker_reason=checked[1],
                                                 source_accepted=bool(accepted)))
            if accepted:
                after = [orient(e) for e in source.elements()]
                checks = check_front_invariants(domain, orient(source.front()), after)
                failed = [k for k, v in checks.items() if k != 'ok' and not v]
                added = after[len(elements):]
                if candidate is not None and (not added or added[0].shape != candidate.shape or
                                              not any(np.allclose(added[0], np.roll(candidate, shift, axis=0),
                                                                  rtol=0, atol=tolerance)
                                                      for shift in range(len(candidate)))):
                    failed.append('candidate_not_added')
                if len(added) > 1:
                    report.closing.append(dict(record, extra_elements=[e.tolist() for e in added[1:]]))
                if failed:
                    report.counts['invariant_failures'] += 1
                    report.invariant_failures.append(dict(record, failed_keys=failed))
            if done:
                break
    return report
