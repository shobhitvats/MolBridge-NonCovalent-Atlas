"""Tests for performance.timing utilities."""
import time
from performance.timing import time_block, time_function, TIMINGS

def test_time_block_and_snapshot():
    TIMINGS.clear()
    with time_block("dummy", items=5):
        time.sleep(0.01)
    snap = TIMINGS.snapshot()
    assert "dummy" in snap
    assert snap["dummy"]["calls"] == 1
    assert snap["dummy"]["total_items"] == 5
    assert snap["dummy"]["total_time"] > 0


def test_time_function_decorator():
    TIMINGS.clear()
    @time_function(items_attr="n")
    def produce():
        class R:
            n = 7
        time.sleep(0.005)
        return R()
    produce()
    snap = TIMINGS.snapshot()
    assert "produce" in snap
    assert snap["produce"]["calls"] == 1
    assert snap["produce"]["total_items"] == 7
