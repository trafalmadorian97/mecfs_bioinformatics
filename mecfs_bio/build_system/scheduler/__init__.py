"""
Given target assets, the scheduler is responsible for traversing the dependency graph to reach those target assets.

Delegates the actual work of materializing an individual asset up to date to the rebuilder.
"""
