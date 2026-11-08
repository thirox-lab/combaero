import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import TypedDict


class ExampleResult(TypedDict):
    name: str
    success: bool
    duration: float
    error: str | None


# Examples that are known-broken for a reason that is not a quick fix, and why.
# These are reported as SKIP rather than failing the run, so that CI stays a
# real signal instead of a permanently red check nobody reads. An entry here is
# tracked debt, not a dismissal: each one names the issue that owns it, and the
# entry comes out the moment the example runs again.
KNOWN_BROKEN: dict[str, str] = {}

# The default per-example limit [s], and the examples that need longer, with
# why. network_from_grid.py solves three growing coupled grids (its last, 4 x 4
# with walls, is the heaviest network any example builds); main ran it in
# 82-116 s on CI against the 120 s default and timed out once, and the exact
# wall composition relay (#717) adds about 10% more relay passes.
DEFAULT_TIMEOUT = 120
SLOW_EXAMPLES: dict[str, int] = {"network_from_grid.py": 300, "hylton1983_c3x_vane.py": 600}


def notebook_source(path: Path) -> str:
    """The code cells of a notebook, joined into one script."""
    cells = json.loads(path.read_text())["cells"]
    return "\n\n".join("".join(c["source"]) for c in cells if c["cell_type"] == "code")


def run_example(example_path: Path, timeout: int = DEFAULT_TIMEOUT) -> ExampleResult:
    """Run a single example script and return the result."""
    start_time = time.perf_counter()
    try:
        # Run with the same python interpreter and headless plotting
        env = os.environ.copy()
        env["COMBAERO_SAVE_PLOTS"] = "1"

        if example_path.suffix == ".ipynb":
            # A notebook runs as its code cells joined into one script, from
            # the repository root, with no Jupyter dependency. Cell outputs
            # are not checked; a cell that raises fails the run. Agg keeps
            # plt.show() from opening a window.
            env["MPLBACKEND"] = "Agg"
            cmd = [sys.executable, "-c", notebook_source(example_path)]
        else:
            cmd = [sys.executable, str(example_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
        duration = time.perf_counter() - start_time
        if result.returncode == 0:
            return {
                "name": example_path.name,
                "success": True,
                "duration": duration,
                "error": None,
            }
        else:
            return {
                "name": example_path.name,
                "success": False,
                "duration": duration,
                "error": result.stderr.strip() or "Process exited with non-zero return code.",
            }
    except subprocess.TimeoutExpired:
        duration = time.perf_counter() - start_time
        return {
            "name": example_path.name,
            "success": False,
            "duration": duration,
            "error": f"Timeout expired after {timeout}s",
        }
    except Exception as e:
        duration = time.perf_counter() - start_time
        return {
            "name": example_path.name,
            "success": False,
            "duration": duration,
            "error": str(e),
        }


def main() -> None:
    """Discover and run all Python examples."""
    examples_dir = Path("python/examples")
    if not examples_dir.exists():
        print(f"Error: {examples_dir} not found.")
        sys.exit(1)

    # Discover all .py files, excluding utilities, the notebooks (*.ipynb),
    # and the worked examples from papers (#793) in papers/
    example_files = sorted([f for f in examples_dir.glob("*.py") if f.name != "plot_utils.py"])
    example_files += sorted(examples_dir.glob("*.ipynb"))
    example_files += sorted((examples_dir / "papers").glob("*.py"))

    print(f"Found {len(example_files)} examples in {examples_dir}")
    print("-" * 60)

    results: list[ExampleResult] = []

    for i, example_file in enumerate(example_files, 1):
        print(f"[{i}/{len(example_files)}] Running {example_file.name}...", end=" ", flush=True)
        res = run_example(example_file, SLOW_EXAMPLES.get(example_file.name, DEFAULT_TIMEOUT))
        results.append(res)

        if res["name"] in KNOWN_BROKEN:
            # Known-broken examples are still run, so that one which starts
            # working again is noticed instead of sitting on the list forever.
            print("PASSED (was expected to fail)" if res["success"] else "SKIP (known broken)")
        elif res["success"]:
            print(f"SUCCESS ({res['duration']:.2f}s)")
        else:
            print("FAILED")

    print("\n" + "=" * 60)
    print("SUMMARY REPORT")
    print("=" * 60)

    successful = [r for r in results if r["success"] and r["name"] not in KNOWN_BROKEN]
    failed = [r for r in results if not r["success"] and r["name"] not in KNOWN_BROKEN]
    skipped = [r for r in results if not r["success"] and r["name"] in KNOWN_BROKEN]
    unexpected_pass = [r for r in results if r["success"] and r["name"] in KNOWN_BROKEN]

    total_time = sum(r["duration"] for r in results)

    print(f"Total Examples: {len(results)}")
    print(f"Successful:     {len(successful)}")
    print(f"Failed:         {len(failed)}")
    print(f"Skipped:        {len(skipped)} (known broken)")
    print(f"Total Wall Time: {total_time:.2f}s")
    print("-" * 60)

    if skipped:
        print("\nSKIPPED (known broken, tracked):")
        for res in skipped:
            print(f"  {res['name']}")
            print(f"    {KNOWN_BROKEN[res['name']]}")
        print("-" * 60)

    if failed:
        print("\nFAILED EXAMPLES DIAGNOSIS:")
        for res in failed:
            print(f"\n--- {res['name']} ---")
            print(f"Error:\n{res['error']}")
        print("-" * 60)

    if unexpected_pass:
        print("\nThese examples are listed in KNOWN_BROKEN but now pass.")
        print("Remove their entries from scripts/run_examples.py:")
        for res in unexpected_pass:
            print(f"  {res['name']}")
        print("-" * 60)

    if failed or unexpected_pass:
        sys.exit(1)

    print("\nAll examples ran successfully!")
    sys.exit(0)


if __name__ == "__main__":
    main()
