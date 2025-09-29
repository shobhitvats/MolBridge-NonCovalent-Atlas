from analysis.keys import resolve_key
from utils.normalization import normalize_detector_output


def test_resolve_key_aliases():
    assert resolve_key('hydrogenbond') == 'hydrogen_bond'
    assert resolve_key('HYDROGENBOND') == 'hydrogen_bond'
    assert resolve_key('ch_pi') == 'ch_pi'
    assert resolve_key('npistar') == 'n_pi_star'


def test_normalize_detector_output_list_mixed():
    raw = [{'distance': 3.1}, object()]  # second will be coerced
    norm = normalize_detector_output('hydrogen_bond', raw)
    assert len(norm) == 2
    assert all('type' in n for n in norm)


def test_normalize_detector_output_single():
    norm = normalize_detector_output('pi_pi', {'angle': 25})
    assert norm[0]['type'] == 'pi_pi'
