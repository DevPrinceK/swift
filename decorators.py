
"""
Decorator to calculate the time a function takes to execute.
"""

from time import perf_counter


def timer(func):
    """
    Decorator to measure the execution time of a function.
    """
    def wrapper(*args, **kwargs):
        start_time = perf_counter()
        result = func(*args, **kwargs)
        end_time = perf_counter()
        print(f"Execution time: {end_time - start_time:.4f} seconds")
        return result
    return wrapper