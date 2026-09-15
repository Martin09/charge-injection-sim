"""Solve and export the Figure 4(b) benchmark or a saved input snapshot."""

import argparse
from datetime import UTC, datetime
from pathlib import Path

from charge_injection_sim import load_inputs, save_run, solve_all_cases


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="benchmark TOML or saved literature-unit inputs.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="new output directory (default: outputs/figure4b-<UTC timestamp>)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = load_inputs(args.config)
    resolved = inputs.to_si()
    results = solve_all_cases(resolved)
    output = args.output or Path("outputs") / datetime.now(UTC).strftime("figure4b-%Y%m%dT%H%M%SZ")
    saved_path = save_run(output, inputs, resolved, results)
    print(saved_path)


if __name__ == "__main__":
    main()
