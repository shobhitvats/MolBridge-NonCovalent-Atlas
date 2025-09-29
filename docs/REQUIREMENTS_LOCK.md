Requirements Lock Generation
============================

Purpose
-------
Provide a deterministic, reproducible dependency set for performance regression tracking and deployment.

Tools
-----
You may use either `pip-tools` (preferred) or `uv` (fast resolver) depending on environment constraints.

Option A: pip-tools
-------------------
1. Install pip-tools (one time):
   pip install pip-tools
2. Maintain abstract deps only in `pyproject.toml` (project section) and minimal `requirements.in` (if needed for extras/testing).
3. Compile lock:
   pip-compile --resolver=backtracking --generate-hashes -o requirements.lock.txt pyproject.toml
4. For dev/test extras:
   pip-compile --extra test --generate-hashes -o requirements.test.lock.txt pyproject.toml
5. Sync environment (optional CI step):
   pip-sync requirements.lock.txt requirements.test.lock.txt

Option B: uv
------------
1. Install:
   curl -LsSf https://astral.sh/uv/install.sh | sh
2. Export lock:
   uv pip compile pyproject.toml -o requirements.lock.txt
3. Install from lock:
   uv pip install -r requirements.lock.txt

CI Recommendations
------------------
Add a workflow step that fails if `pip-compile` would modify the existing lock (detect drift):
  - Run compile in CI
  - git diff --exit-code requirements.lock.txt

Performance Baseline Coupling
-----------------------------
Tie a specific lock commit hash to the golden dataset manifest to guarantee geometry & numeric determinism for regression comparisons.

Update Policy
-------------
- Routine: every 4 weeks or when adding new detector dependencies.
- Emergency: security patch releases (CVE) as required.

Checksum Verification
---------------------
Use `--generate-hashes` so supply chain integrity is enforced during install (`pip install --require-hashes -r requirements.lock.txt`).

Regeneration Script (Optional)
------------------------------
Create `scripts/update_lock.sh`:
  #!/usr/bin/env bash
  set -euo pipefail
  pip install --upgrade pip pip-tools
  pip-compile --resolver=backtracking --generate-hashes -o requirements.lock.txt pyproject.toml
  pip-compile --extra test --generate-hashes -o requirements.test.lock.txt pyproject.toml
  echo "Locks updated. Review and commit." 

Future Enhancements
-------------------
- Add SBOM generation (cyclonedx-py).
- Integrate vulnerability scanning (pip-audit) against locked set.
