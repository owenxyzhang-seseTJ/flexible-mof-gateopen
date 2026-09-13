#!/usr/bin/env python3
"""Download a public CIF with optional SHA-256 provenance checking."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from urllib.request import Request, urlopen


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sha256", default=None)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    request = Request(args.url, headers={"User-Agent": "flexible-mof-gateopen/0.1"})
    with urlopen(request, timeout=60) as response:
        args.output.write_bytes(response.read())
    digest = sha256(args.output)
    if args.sha256 and digest.lower() != args.sha256.lower():
        raise SystemExit(f"SHA-256 mismatch: expected {args.sha256}, got {digest}")
    print(f"downloaded {args.output}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
