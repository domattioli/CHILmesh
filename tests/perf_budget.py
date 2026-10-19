"""Wall-clock budgets for timing tests.

Shared CI runners are noisy, so every timing assert widens its limit 3x
when ``CI`` is set. Local runs keep the original limit.
"""
import os

CI_SLACK = 3.0 if os.environ.get("CI") else 1.0


def budget(seconds):
    """Return the time limit for this environment."""
    return seconds * CI_SLACK
