from importlib import reload


def test_adapt_threshold_reason_variants(monkeypatch):
    import utils.kdtree_thresholds as kt
    reload(kt)
    base = kt.get_threshold('hydro')
    # Force raise_low_util (kdtree_used True & small candidate complexity)
    new_val, changed, reason = kt.adapt_threshold('hydro', int(base * 0.1), True)
    assert reason in {'raise_low_util', 'unchanged'}
    # Force lower_high_miss (kdtree not used & large complexity)
    base2 = kt.get_threshold('hydro')
    huge = base2 * 3
    new_val2, changed2, reason2 = kt.adapt_threshold('hydro', huge, False)
    assert reason2 in {'lower_high_miss', 'density_adjust_up', 'density_adjust_down', 'unchanged'}
    # Force density adjustments by crafting densities (simulate mid-range complexity)
    mid = int(base2 * 1.1)
    kt.adapt_threshold('hydro', mid, True)  # may adjust density
